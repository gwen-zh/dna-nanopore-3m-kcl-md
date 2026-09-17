"""Bounded read-only dA pore analysis and coordinate-derived movie inputs."""
import argparse
import hashlib
import itertools
import json
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from trajectory import psf,edges,whole,binary,DCD,write_dcd,write_psf,write_pdb

AMINO={'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','HSD','HSE','HSP','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}


def sha(path,start=0,limit=None):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        f.seek(start)
        left=(Path(path).stat().st_size-start) if limit is None else limit
        while left:
            block=f.read(min(1048576,left))
            if not block: raise EOFError('Hash input shortened')
            h.update(block);left-=len(block)
    return h.hexdigest()


def rotation(x,target):
    center=x.mean(0);reference=target.mean(0)
    u,_,vt=np.linalg.svd((x-center).T@(target-reference))
    r=u@np.diag([1,1,np.linalg.det(u@vt)])@vt
    return center,r,reference


def run(args):
    root=args.workspace.resolve();pore=args.pore;out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    job={'7ahl':'57323201','3b07':'57323197'}[pore]
    source=root/f'pore-fix-20260913/polyda40-3m-kcl-{pore}/run-{job}'
    system=root/f'polyda40-3m-kcl-{pore}/formal-50ns'
    events=[json.loads(s) for s in (source/'events.jsonl').read_text().splitlines()]
    completed=any(e['event']=='workflow_complete' for e in events)
    atoms,bonds=psf(system/'system.psf');n=len(atoms)
    dna=np.flatnonzero(atoms.segid.to_numpy()=='AN1')
    protein=np.flatnonzero(atoms.resname.isin(AMINO).to_numpy())
    heavy=atoms.mass.to_numpy()>2
    ph=protein[heavy[protein]];dh=dna[heavy[dna]]
    ca=protein[atoms.iloc[protein].name.to_numpy()=='CA']
    lipidp=np.flatnonzero((atoms.resname.to_numpy()=='POPC')&(atoms.name.to_numpy()=='P'))
    if len(dna)!=1279 or set(atoms.iloc[dna].resname)!={'ADE'}: raise ValueError('Expected poly(dA)40')
    mapping=np.full(n,-1);mapping[dna]=np.arange(len(dna))
    db=mapping[bonds[np.isin(bonds,dna).all(1)]];dtree=edges(len(dna),db)
    protein_trees=[]
    for seg in atoms.iloc[protein].segid.unique():
        ids=protein[atoms.iloc[protein].segid.to_numpy()==seg]
        inv=np.full(n,-1);inv[ids]=np.arange(len(ids))
        protein_trees.append((ids,edges(len(ids),inv[bonds[np.isin(bonds,ids).all(1)]])))
    reference=binary(source/'eq3.coor',n)
    refbox=np.loadtxt(source/'eq3.xsc',comments='#')[[1,5,9]]
    for ids,tree in protein_trees: reference[ids]=whole(reference[ids],refbox,tree)
    ca_reference=reference[ca].copy()
    cis_z=float(ca_reference[:,2].max());trans_z=float(ca_reference[:,2].min())
    pore_xy=ca_reference[:,:2].mean(0)
    c1=np.flatnonzero(atoms.iloc[dna].name.to_numpy()=="C1'")
    lookup={(int(a.resid),a['name']):i for i,a in atoms.iloc[dna].reset_index(drop=True).iterrows()}
    rings=np.array([[lookup[r,k] for k in ['N9','C4','C8']] for r in range(1,41)])
    bases=np.array([[lookup[r,k] for k in ['N9','C8','N7','C5','C6','N1','C2','N3','C4']] for r in range(1,41)])
    selected=np.sort(np.r_[ph,dna,lipidp]);inv=np.full(n,-1);inv[selected]=np.arange(len(selected))
    visual_bonds=inv[bonds[np.isin(bonds,selected).all(1)]]
    dcd=DCD(source/'prod5.dcd',n,completed)
    if dcd.first!=25000 or dcd.stride!=25000 or (completed and dcd.count!=200): raise ValueError('Unexpected production time axis')
    last_center=reference[dna].mean(0)
    rows=[];positions=[];contacts_by_residue=[];movie=[];dna_coords=[]
    occupancy={};image_min=float('inf');worst_bond=0.

    def consume(step,raw,box):
        nonlocal last_center,image_min,worst_bond
        x=whole(raw[dna],box,dtree)
        x+=box*np.rint((last_center-x.mean(0))/box);last_center=x.mean(0)
        lengths=np.linalg.norm(x[db[:,0]]-x[db[:,1]],axis=1)
        if lengths.min()<=.5 or lengths.max()>=2.2: raise ValueError('Broken complete DNA')
        worst_bond=max(worst_bond,float(lengths.max()))
        protein_xyz=raw.copy()
        for ids,tree in protein_trees:
            part=whole(raw[ids],box,tree)
            part+=box*np.rint((reference[ids].mean(0)-part.mean(0))/box)
            protein_xyz[ids]=part
        center,r,refcenter=rotation(protein_xyz[ca],ca_reference)
        aligned=(x-center)@r+refcenter
        mass=atoms.mass.to_numpy()[dh];xh=aligned[heavy[dna]]
        com=np.average(xh,axis=0,weights=mass)
        rg=float(np.sqrt(np.average(np.sum((xh-com)**2,axis=1),weights=mass)))
        tree=cKDTree(np.mod(raw[ph],box),boxsize=box)
        nearest=float(tree.query(np.mod(x[heavy[dna]],box))[0].min())
        neighbours=tree.query_ball_point(np.mod(x[heavy[dna]],box),3.5)
        dcontact=np.array([bool(v) for v in neighbours])
        pcontact=np.unique([j for v in neighbours for j in v]).astype(int)
        dres=np.unique(atoms.resid.to_numpy()[dh[dcontact]])
        contact_row=np.isin(np.arange(1,41),dres);contacts_by_residue.append(contact_row)
        for _,a in atoms.iloc[ph[pcontact]].drop_duplicates(['segid','resid']).iterrows():
            key=(a.segid,int(a.resid),a.resname);occupancy[key]=occupancy.get(key,0)+(step>0)
        heavy_whole=x[heavy[dna]]
        if np.any(np.ptp(heavy_whole,axis=0)>=box): raise ValueError('DNA exceeds image-audit assumption')
        itree=cKDTree(heavy_whole)
        gap=min(float(itree.query(heavy_whole+np.array(s)*box)[0].min())
                for s in itertools.product([-1,0,1],repeat=3) if s!=(0,0,0))
        image_min=min(image_min,gap)
        bc=aligned[bases].mean(1);rp=aligned[rings]
        normals=np.cross(rp[:,1]-rp[:,0],rp[:,2]-rp[:,0]);normals/=np.linalg.norm(normals,axis=1)[:,None]
        v=bc[1:]-bc[:-1];distance=np.linalg.norm(v,axis=1)
        a=np.abs(np.sum(normals[:-1]*normals[1:],axis=1));hi=np.abs(np.sum(v*normals[:-1],axis=1));hj=np.abs(np.sum(v*normals[1:],axis=1))
        stack=(distance<=5.5)&(a>=np.cos(np.pi/4))&(hi>=2)&(hj>=2)&(hi<=4.5)&(hj<=4.5)&(np.sqrt(np.maximum(0,distance**2-hi**2))<=3.5)&(np.sqrt(np.maximum(0,distance**2-hj**2))<=3.5)
        lead=aligned[c1[0]];tail=aligned[c1[-1]]
        below=int((aligned[c1,2]<trans_z).sum())
        rows.append(dict(step=step,time_ns=step*1e-6,rg_A=rg,end_to_end_A=float(np.linalg.norm(tail-lead)),
                         lead_z_A=float(lead[2]),tail_z_A=float(tail[2]),dna_COM_z_A=float(com[2]),
                         lead_radius_A=float(np.linalg.norm(lead[:2]-pore_xy)),lead_above_cis_CA_plane_A=float(lead[2]-cis_z),
                         c1_below_trans_CA_plane=below,minimum_DNA_protein_heavy_A=nearest,
                         contact_pairs_3p5A=sum(map(len,neighbours)),contacting_DNA_residues_3p5A=len(dres),
                         contacting_protein_atoms_3p5A=len(pcontact),stacked_adjacent_pairs_ring_centers=int(stack.sum()),
                         dna_max_bond_A=float(lengths.max()),nearest_DNA_image_A=gap))
        display=protein_xyz[selected].copy()
        # Put lipid headgroups into the pore-centred image, without creating
        # membrane-wide bonds. Only DNA uses temporal continuity.
        lpidx=inv[lipidp]
        display[lpidx]=raw[lipidp]+box*np.rint((refcenter-raw[lipidp])/box)
        display=(display-center)@r+refcenter;display[inv[dna]]=aligned
        movie.append(display.astype(np.float32));dna_coords.append(aligned.astype(np.float32));positions.append(aligned[c1,2])
        if len(rows)%25==0: print(pore,'processed',len(rows),'frames',flush=True)

    consume(0,reference.copy(),refbox)
    for step,raw,box in dcd: consume(step,raw,box)
    table=pd.DataFrame(rows);table.to_csv(out/'timeseries.csv',index=False)
    z=np.array(positions);cm=np.array(contacts_by_residue)
    pd.DataFrame(z,columns=[f'nt_{i}' for i in range(1,41)]).assign(time_ns=table.time_ns).to_csv(out/'nucleotide_z.csv',index=False)
    pd.DataFrame(cm.astype(int),columns=[f'nt_{i}' for i in range(1,41)]).assign(time_ns=table.time_ns).to_csv(out/'nucleotide_contacts.csv',index=False)
    records=[dict(chain=k[0],resid=k[1],resname=k[2],contact_frame_fraction=v/dcd.count,contact_frames=v) for k,v in occupancy.items() if v]
    pd.DataFrame(records,columns=['chain','resid','resname','contact_frame_fraction','contact_frames']).sort_values('contact_frame_fraction',ascending=False).to_csv(out/'protein_contact_occupancy.csv',index=False)
    np.savez_compressed(out/'dna_coordinates.npz',xyz=np.array(dna_coords),time_ns=table.time_ns.to_numpy())
    vcoords=np.array(movie);v_atoms=atoms.iloc[selected].reset_index(drop=True)
    write_psf(out/'visualization.psf',v_atoms,visual_bonds)
    write_pdb(out/'start.pdb',v_atoms,vcoords[0]);write_pdb(out/'endpoint.pdb',v_atoms,vcoords[-1])
    write_dcd(out/'visualization.dcd',vcoords,dcd.header)
    write_pdb(out/'dna_start.pdb',atoms.iloc[dna].reset_index(drop=True),dna_coords[0])
    write_pdb(out/'dna_endpoint.pdb',atoms.iloc[dna].reset_index(drop=True),dna_coords[-1])
    final=rows[-1];first=rows[0]
    summary=dict(pore=pore,job_id=job,status='complete' if completed else 'preview_in_progress',
                 captured_utc=datetime.now(timezone.utc).isoformat(),production_ns=final['time_ns'],target_ns=5,
                 production_frames=dcd.count,includes_time_zero_checkpoint=True,frame_spacing_ps=25,timestep_fs=1,
                 protein_alignment='Rigid proper fit of all protein CA atoms to production-start checkpoint; shared transform for all displayed atoms',
                 cis_CA_plane_z_A=cis_z,trans_CA_plane_z_A=trans_z,initial=first,endpoint=final,
                 lead_axial_advance_toward_minus_z_A=first['lead_z_A']-final['lead_z_A'],
                 all40_c1_below_trans_plane_observed=bool(np.any(table.c1_below_trans_CA_plane==40)),
                 minimum_DNA_image_A=image_min,maximum_DNA_bond_A=worst_bond,
                 average_production_contact_pairs_3p5A=float(table.contact_pairs_3p5A.iloc[1:].mean()),
                 average_production_contacting_DNA_residues=float(table.contacting_DNA_residues_3p5A.iloc[1:].mean()),
                 source_atom_count=n,visualization_atom_count=len(selected),
                 forcefield_and_protocol='See current run protocol; nominal 3 M preparation, ~1 V box-spanning field; not measured local voltage or salt density',
                 interpretation='Single trajectory. Axial progress and contact counts are geometric observables, not binding energies, ionic current, a rate or proof of equilibrium. Complete 5 ns does not imply full translocation.',
                 source_files={'psf':str((system/'system.psf').relative_to(root)),'dcd':str((source/'prod5.dcd').relative_to(root)),
                               'time_zero':str((source/'eq3.coor').relative_to(root))},
                 source_psf_sha256=sha(system/'system.psf'),source_DCD_header_bytes=dcd.offset,
                 source_DCD_snapshot_byte_limit=dcd.byte_limit,
                 source_DCD_frame_payload_sha256=sha(source/'prod5.dcd',dcd.offset,dcd.count*dcd.framebytes),
                 source_DCD_full_sha256=sha(source/'prod5.dcd') if completed else None)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (out/'render.json').write_text(json.dumps(dict(pore=pore,center=[float(pore_xy[0]),float(pore_xy[1]),float((vcoords[:,:,2].max()+vcoords[:,:,2].min())/2)],
              z_bounds=[float(vcoords[:,:,2].min()),float(vcoords[:,:,2].max())],frames=len(vcoords),status=summary['status']),indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(11,7),constrained_layout=True)
    t=table.time_ns
    axes[0,0].plot(t,table.lead_z_A,label="5' C1'",color='#eb8b22');axes[0,0].plot(t,table.tail_z_A,label="3' C1'",color='#3587ae')
    axes[0,0].axhline(cis_z,ls='--',color='#999',label='Protein CA cis plane');axes[0,0].axhline(trans_z,ls=':',color='#999',label='Protein CA trans plane')
    axes[0,0].set(ylabel='Pore-aligned z (A)',title='Axial position (not a lumen-occupancy test)');axes[0,0].legend(fontsize=8)
    axes[0,1].plot(t,table.rg_A,color='#d97919');axes[0,1].set(ylabel='Heavy-atom Rg (A)',title='DNA conformation')
    axes[1,0].plot(t,table.contacting_DNA_residues_3p5A,color='#9256a5');axes[1,0].set(ylabel='Contacting nucleotides',title='Periodic DNA-protein heavy contacts <3.5 A',ylim=(-.5,40.5))
    im=axes[1,1].imshow(z.T,origin='lower',aspect='auto',extent=[0,t.iloc[-1],.5,40.5],cmap='viridis');fig.colorbar(im,ax=axes[1,1],label='C1\' z (A)')
    axes[1,1].set(ylabel='Nucleotide',title='Axial coordinates of all 40 nucleotides')
    for ax in axes.flat: ax.set_xlabel('Production time (ns)')
    fig.suptitle(f'poly(dA)40 / {pore.upper()} | {t.iloc[-1]:.3f} ns | '+('completed segment' if completed else 'in-progress preview'))
    fig.savefig(out/'overview.png',dpi=160);fig.savefig(out/'overview.pdf');plt.close(fig)
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workspace',type=Path,required=True);p.add_argument('--pore',choices=['7ahl','3b07'],required=True);p.add_argument('--out',type=Path,required=True)
    run(p.parse_args())
