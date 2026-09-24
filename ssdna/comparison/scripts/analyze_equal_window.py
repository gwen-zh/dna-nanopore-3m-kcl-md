#!/usr/bin/env python3
"""Matched-age 5 ns structural analysis for dA40 and dT40.

Both sequences use the 0.01–5.00 ns production window after the same heating
and equilibration schedule. Coordinates are reconstructed through PSF bonds
before any molecular geometry is measured.
"""
import argparse
import collections
import csv
import hashlib
import json
import math
import struct
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist


SUGAR_BACKBONE = {
    'P','O1P','O2P',"O5'","C5'","H5'","H5''","C4'","H4'","O4'",
    "C1'","H1'","C2'","H2'","H2''","C3'","H3'","O3'",'H5T','H3T'
}
COLORS={'dA40':'#c23b22','dT40':'#2468a2'}


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(4*1024*1024),b''):h.update(block)
    return h.hexdigest()


def read_psf(path):
    lines=Path(path).read_text().splitlines()
    start=next(i for i,line in enumerate(lines) if '!NATOM' in line)
    n=int(lines[start].split()[0]);rows=[line.split() for line in lines[start+1:start+1+n]]
    atoms=pd.DataFrame(rows,columns=['id','segid','resid','resname','name','type','charge','mass','unused'])
    atoms[['id','resid']]=atoms[['id','resid']].astype(int)
    atoms[['charge','mass']]=atoms[['charge','mass']].astype(float)
    start=next(i for i,line in enumerate(lines) if '!NBOND' in line);nb=int(lines[start].split()[0]);values=[]
    for line in lines[start+1:]:
        values.extend(map(int,line.split()))
        if len(values)>=2*nb:break
    return atoms,np.asarray(values,dtype=int).reshape(-1,2)-1


class DCD:
    def _record(self):
        marker=self.stream.read(4)
        if len(marker)!=4:raise EOFError
        size=struct.unpack('<i',marker)[0];payload=self.stream.read(size)
        if len(payload)!=size or self.stream.read(4)!=marker:raise EOFError('truncated DCD record')
        return payload

    def __init__(self,path):
        self.path=Path(path);self.stream=self.path.open('rb');header=self._record()
        if len(header)!=84 or header[:4]!=b'CORD':raise ValueError('unsupported DCD')
        control=np.frombuffer(header[4:],dtype='<i4')
        self.nframes,self.firststep,self.stride=map(int,control[:3])
        if control[8]!=0 or control[10]!=1:raise ValueError('fixed atoms or missing cell unsupported')
        self._record();self.natoms=struct.unpack('<i',self._record())[0];self.start=self.stream.tell()
        self.framebytes=56+3*(8+4*self.natoms)
        expected=self.start+self.nframes*self.framebytes
        if self.path.stat().st_size!=expected:raise ValueError('DCD size/header mismatch')

    def frame(self,index):
        self.stream.seek(self.start+index*self.framebytes);cell=np.frombuffer(self._record(),dtype='<f8')
        xyz=np.column_stack([np.frombuffer(self._record(),dtype='<f4') for _ in range(3)]).astype(float)
        return xyz,cell[[0,2,5]].copy()

    def close(self):self.stream.close()


def minimum_image(delta,box):return delta-box*np.rint(delta/box)


def traversal(n,bonds):
    adjacency=[[] for _ in range(n)]
    for i,j in bonds:adjacency[i].append(j);adjacency[j].append(i)
    seen={0};queue=collections.deque([0]);edges=[]
    while queue:
        i=queue.popleft()
        for j in adjacency[i]:
            if j not in seen:seen.add(j);queue.append(j);edges.append((i,j))
    if len(seen)!=n:raise ValueError('DNA is not one connected fragment')
    return edges


def join(raw,box,edges):
    result=raw.copy()
    for i,j in edges:result[j]=result[i]+minimum_image(raw[j]-raw[i],box)
    return result


def unit(value):
    norm=np.linalg.norm(value,axis=-1,keepdims=True)
    return value/np.maximum(norm,1e-12)


def dihedral(points):
    a,b,c,d=np.moveaxis(points,-2,0);axis=unit(c-b)
    v=a-b;w=d-c;v-=np.sum(v*axis,axis=-1,keepdims=True)*axis;w-=np.sum(w*axis,axis=-1,keepdims=True)*axis
    return np.degrees(np.arctan2(np.sum(np.cross(axis,v)*w,axis=-1),np.sum(v*w,axis=-1)))


def kabsch(x,target):
    xx=x-x.mean(0);yy=target-target.mean(0);u,_,vt=np.linalg.svd(xx.T@yy)
    rotation=u@np.diag([1,1,np.sign(np.linalg.det(u@vt))])@vt
    return xx@rotation+target.mean(0)


def circular_mean_deg(values,axis=0):
    rad=np.deg2rad(values)
    return np.rad2deg(np.arctan2(np.nanmean(np.sin(rad),axis=axis),np.nanmean(np.cos(rad),axis=axis)))


def source_setup(label,psf,dcd_path,frames,start_frame):
    atoms,bonds=read_psf(psf)
    dna_ids=atoms.index[atoms.resname.isin(['ADE','THY'])].to_numpy()
    if len(dna_ids)!=1279 or not np.array_equal(dna_ids,np.arange(1279)):
        raise ValueError(f'{label}: expected first 1279 atoms to be dA40/dT40')
    dna=atoms.iloc[dna_ids].reset_index(drop=True)
    if set(dna.resid)!=set(range(1,41)) or len(set(dna.resname))!=1:raise ValueError(label+': invalid DNA residues')
    dbonds=bonds[np.all(np.isin(bonds,dna_ids),axis=1)]
    edges=traversal(len(dna),dbonds)
    dcd=DCD(dcd_path)
    if dcd.natoms!=len(atoms) or start_frame<0 or start_frame+frames>dcd.nframes:
        raise ValueError(label+': incompatible DCD or frame window')
    selected=np.arange(start_frame,start_frame+frames)
    return dict(label=label,psf=Path(psf),dcd=dcd,dcd_path=Path(dcd_path),atoms=atoms,dna=dna,
                dna_ids=dna_ids,dbonds=dbonds,edges=edges,selected=selected)


def analyze_one(info,out,environment_stride=10):
    label=info['label'];atoms=info['atoms'];dna=info['dna'];dcd=info['dcd'];selected=info['selected']
    mass=dna.mass.to_numpy();heavy=np.flatnonzero(mass>2);heavy_mass=mass[heavy]
    residue=dna.resid.to_numpy()-1;lookup={(int(row.resid),row['name']):i for i,row in dna.iterrows()}
    adenine=dna.resname.iloc[0]=='ADE';sequence='poly(dA)40' if adenine else 'poly(dT)40'
    base_ids=[np.flatnonzero((residue==r)&~dna.name.isin(SUGAR_BACKBONE)) for r in range(40)]
    if not all(len(ids)==14 for ids in base_ids):raise ValueError(label+': unexpected base atom count')
    base_heavy=[np.array([i for i in ids if mass[i]>2]) for ids in base_ids]
    c1=np.array([lookup[r,"C1'"] for r in range(1,41)])
    plane_names=['N9','C4','C8'] if adenine else ['N1','C2','C6']
    plane=np.array([[lookup[r,name] for name in plane_names] for r in range(1,41)])
    chi_names=["O4'","C1'",'N9','C4'] if adenine else ["O4'","C1'",'N1','C2']
    chi_ids=np.array([[lookup[r,name] for name in chi_names] for r in range(1,41)])
    donors=[]
    for r in range(1,41):
        names=[('N6','H61'),('N6','H62')] if adenine else [('N3','H3')]
        donors.extend((lookup[r,a],lookup[r,h]) for a,h in names)
    donors=np.asarray(donors)
    accept_names=(['N1','N3','N7'] if adenine else ['O2','O4'])+['O1P','O2P',"O3'","O5'","O4'"]
    acceptors=np.flatnonzero(dna.name.isin(accept_names));base_flat=np.concatenate(base_heavy)
    donor_base=np.isin(donors[:,0],base_flat);acceptor_base=np.isin(acceptors,base_flat)
    donor_acceptor_mask=residue[donors[:,0],None]!=residue[acceptors][None,:]
    k_ids=atoms.index[atoms.name=='POT'].to_numpy();cl_ids=atoms.index[atoms.name=='CLA'].to_numpy()
    water_ids=atoms.index[(atoms.resname=='TIP3')&(atoms.name=='OH2')].to_numpy()
    rows=[];local=[];environment=[];contact_maps=[];stack_maps=[];hbond_maps=[];coordinates=[];base_centers=[]
    maximum_bond=0.;source_times=[]
    for count,index in enumerate(selected,1):
        raw,box=dcd.frame(int(index));x=join(raw[info['dna_ids']],box,info['edges'])
        bond=np.linalg.norm(x[info['dbonds'][:,0]]-x[info['dbonds'][:,1]],axis=1)
        if bond.min()<=.5 or bond.max()>=2.2:raise ValueError(label+': DNA bond gate failed')
        maximum_bond=max(maximum_bond,float(bond.max()))
        source_time=(dcd.firststep+index*dcd.stride)*2e-6;source_times.append(source_time)
        window_time=count*(dcd.stride*2e-6)
        center=np.average(x[heavy],axis=0,weights=heavy_mass);q=x[heavy]-center
        eigen=np.linalg.eigvalsh((q*heavy_mass[:,None]).T@q/heavy_mass.sum());rg2=eigen.sum()
        c=x[c1];contour=np.linalg.norm(np.diff(c,axis=0),axis=1).sum();end=np.linalg.norm(c[-1]-c[0])
        bc=np.array([np.average(x[ids],axis=0,weights=mass[ids]) for ids in base_ids])
        p=x[plane];normal=unit(np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]))
        delta=bc[None,:,:]-bc[:,None,:];distance=np.linalg.norm(delta,axis=2);alignment=np.abs(normal@normal.T)
        hi=np.abs(np.einsum('ijk,ik->ij',delta,normal));hj=np.abs(np.einsum('ijk,jk->ij',delta,normal))
        li=np.sqrt(np.maximum(0,distance**2-hi**2));lj=np.sqrt(np.maximum(0,distance**2-hj**2))
        strict=(distance<=5.5)&(alignment>=math.cos(math.pi/4))&(hi>=2)&(hj>=2)&(hi<=4.5)&(hj<=4.5)&(li<=3.5)&(lj<=3.5)
        np.fill_diagonal(strict,False);adj=np.diag(strict,1);neighbor=np.r_[False,adj]|np.r_[adj,False]
        longmask=np.abs(np.arange(40)[:,None]-np.arange(40)[None,:])>=3
        chi=dihedral(x[chi_ids]);outward=np.sum(unit(bc-c)*unit(c-center),axis=1)
        tangent=unit(np.vstack([c[1]-c[0],c[2:]-c[:-2],c[-1]-c[-2]]))
        normal_angle=np.degrees(np.arccos(np.clip(np.abs(np.sum(normal*tangent,axis=1)),0,1)))
        contacts=np.zeros((40,40),bool)
        for ai,bi in cKDTree(x[heavy]).query_pairs(4.5):
            ra,rb=residue[heavy[ai]],residue[heavy[bi]]
            if abs(ra-rb)>=3:contacts[ra,rb]=contacts[rb,ra]=True
        dr=x[donors[:,0]];hr=x[donors[:,1]];ar=x[acceptors]
        da=cdist(dr,ar);cosang=np.sum(unit(dr-hr)[:,None,:]*unit(ar[None,:,:]-hr[:,None,:]),axis=-1)
        hb=(da<=3.5)&(cosang<=math.cos(math.radians(150)))&donor_acceptor_mask
        hbbase=hb&donor_base[:,None]&acceptor_base[None,:];hbmap=np.zeros((40,40),bool)
        di,aj=np.where(hbbase);hbmap[residue[donors[di,0]],residue[acceptors[aj]]]=True
        row=dict(sequence=label,window_time_ns=window_time,source_time_ns=source_time,
                 rg_A=math.sqrt(rg2),end_to_end_A=end,c1_contour_A=contour,extension_fraction=end/contour,
                 shape_anisotropy=1.5*np.sum(eigen**2)/rg2**2-.5,
                 adjacent_stacks=int(adj.sum()),nonadjacent_stacks=int(np.triu(strict&longmask,1).sum()),
                 unstacked_neighbor_fraction=float((~neighbor).mean()),unstacked_any_fraction=float((~strict.any(1)).mean()),
                 syn_fraction=float((np.abs(chi)<90).mean()),outward_projection_mean=float(outward.mean()),
                 nonlocal_contact_pairs=int(np.triu(contacts,1).sum()),base_base_hbonds=int(hbbase.sum()),
                 all_interresidue_hbonds=int(hb.sum()),max_DNA_bond_A=float(bond.max()))
        rows.append(row)
        for r in range(40):
            local.append(dict(sequence=label,window_time_ns=window_time,source_time_ns=source_time,resid=r+1,
                              chi_deg=chi[r],syn=int(abs(chi[r])<90),unstacked_neighbors=int(not neighbor[r]),
                              unstacked_any=int(not strict[r].any()),outward_projection=outward[r],
                              base_normal_backbone_deg=normal_angle[r],nonlocal_contact_count=int(contacts[r].sum()),
                              stacked_left=int(r>0 and adj[r-1]),stacked_right=int(r<39 and adj[r])))
        if (count-1)%environment_stride==0 or count==len(selected):
            dna_tree=cKDTree(np.mod(x[heavy],box),boxsize=box)
            envrow=dict(sequence=label,window_time_ns=window_time,source_time_ns=source_time)
            for env_label,ids in [('K',k_ids),('Cl',cl_ids),('waterO',water_ids)]:
                nearest=dna_tree.query(np.mod(raw[ids],box),k=1)[0]
                envrow[f'{env_label}_within_3p5A']=int((nearest<=3.5).sum())
                if env_label!='waterO':envrow[f'{env_label}_box_molar']=len(ids)/(np.prod(box)*6.02214076e-4)
            environment.append(envrow)
        coordinates.append(x.astype(np.float32));base_centers.append(bc.astype(np.float32))
        contact_maps.append(contacts);stack_maps.append(strict);hbond_maps.append(hbmap)
        if count%50==0:print(f'{label}: {count}/{len(selected)} frames',flush=True)
    dcd.close();frame=pd.DataFrame(rows);per=pd.DataFrame(local);coords=np.asarray(coordinates)
    contact_occ=np.asarray(contact_maps).mean(0);stack_occ=np.asarray(stack_maps).mean(0);hb_occ=np.asarray(hbond_maps).mean(0)
    reference=coords[0,c1].astype(float);aligned=np.array([kabsch(v[c1],reference) for v in coords])
    for _ in range(3):reference=np.mean([kabsch(v[c1],reference) for v in coords],axis=0)
    aligned=np.array([kabsch(v[c1],reference) for v in coords]);rmsf=np.sqrt(np.mean(np.sum((aligned-aligned.mean(0))**2,axis=2),axis=0))
    internal=np.array([cdist(v[c1],v[c1]) for v in coords]);target=internal.mean(0)
    representative=int(np.argmin(np.mean((internal-target)**2,axis=(1,2))))
    rep=coords[representative].astype(float);rep-=rep[heavy].mean(0)
    pdb=[]
    for i,atom in dna.iterrows():
        px,py,pz=rep[i];element=next(char for char in atom['name'] if char.isalpha())
        pdb.append(f"ATOM  {i+1:5d} {atom['name']:>4s} {atom.resname:>3s} A{atom.resid:4d}    {px:8.3f}{py:8.3f}{pz:8.3f}{1.:6.2f}{0.:6.2f}      DNA {element:>2s}\n")
    (out/f'{label}_representative.pdb').write_text(f'REMARK {sequence}; real sampled frame; source time {source_times[representative]:.3f} ns\n'+''.join(pdb)+'END\n')
    summary=[]
    local_group=per.groupby('resid')
    for r in range(1,41):
        group=local_group.get_group(r)
        summary.append(dict(sequence=label,resid=r,syn_occupancy=group.syn.mean(),
                            unstacked_neighbor_fraction=group.unstacked_neighbors.mean(),
                            unstacked_any_fraction=group.unstacked_any.mean(),
                            outward_projection_mean=group.outward_projection.mean(),
                            base_normal_backbone_mean_deg=group.base_normal_backbone_deg.mean(),
                            chi_circular_mean_deg=circular_mean_deg(group.chi_deg.to_numpy()),
                            nonlocal_contact_mean=group.nonlocal_contact_count.mean(),C1_RMSF_A=rmsf[r-1]))
    validation=dict(sequence=sequence,label=label,psf=f'{label}_system.psf',dcd=f'{label}_production.dcd',
                    total_DCD_frames=int(info['dcd'].nframes),analyzed_frames=len(selected),
                    selected_frame_range_zero_based=[int(selected[0]),int(selected[-1])],
                    source_time_range_ns=[float(source_times[0]),float(source_times[-1])],
                    analyzed_window_ns=float(frame.window_time_ns.iloc[-1]),maximum_DNA_bond_A=maximum_bond,
                    representative_frame_zero_based=int(selected[representative]),
                    representative_source_time_ns=float(source_times[representative]))
    return dict(frame=frame,local=per,environment=pd.DataFrame(environment),summary=pd.DataFrame(summary),
                contact=contact_occ,stack=stack_occ,hbond=hb_occ,coords=coords,c1=c1,
                base_centers=np.asarray(base_centers),validation=validation)


def save_figure(fig,out,name):
    if not fig.get_constrained_layout():fig.tight_layout()
    fig.savefig(out/(name+'.png'),dpi=220,bbox_inches='tight');fig.savefig(out/(name+'.pdf'),bbox_inches='tight');plt.close(fig)


def plot_all(results,out):
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    metrics=[('rg_A','Rg (Å)'),('end_to_end_A',"5′–3′ C1′ distance (Å)"),('shape_anisotropy','Shape anisotropy'),
             ('adjacent_stacks','Adjacent stacked pairs'),('unstacked_neighbor_fraction','Unstacked-neighbor fraction'),('syn_fraction','syn fraction')]
    fig,axes=plt.subplots(3,2,figsize=(11,10),sharex=True)
    for ax,(key,title) in zip(axes.flat,metrics):
        for label,result in results.items():
            frame=result['frame'];ax.plot(frame.window_time_ns,frame[key],color=COLORS[label],alpha=.18,lw=.7)
            smooth=frame[key].rolling(25,center=True,min_periods=5).mean();ax.plot(frame.window_time_ns,smooth,color=COLORS[label],lw=2,label=label)
        ax.set_ylabel(title);ax.grid(alpha=.2)
    for ax in axes[-1]:ax.set_xlabel('Time in analyzed 5 ns window (ns)')
    axes[0,0].legend(frameon=False);save_figure(fig,out,'01_global_structure')

    profiles=[('syn_occupancy','syn occupancy'),('unstacked_neighbor_fraction','Unstacked-neighbor fraction'),
              ('outward_projection_mean','Outward projection'),('C1_RMSF_A',"C1′ RMSF (Å)")]
    fig,axes=plt.subplots(2,2,figsize=(11,7),sharex=True)
    for ax,(key,title) in zip(axes.flat,profiles):
        for label,result in results.items():ax.plot(result['summary'].resid,result['summary'][key],'-o',ms=2.5,lw=1.4,color=COLORS[label],label=label)
        ax.set_ylabel(title);ax.grid(alpha=.2)
    for ax in axes[-1]:ax.set_xlabel('Residue');axes[0,0].legend(frameon=False)
    save_figure(fig,out,'02_residue_profiles')

    fig,axes=plt.subplots(2,2,figsize=(12,7),sharex=True,sharey=True)
    for col,(label,result) in enumerate(results.items()):
        n=len(result['frame']);syn=result['local'].syn.to_numpy().reshape(n,40).T
        unstack=result['local'].unstacked_neighbors.to_numpy().reshape(n,40).T
        for row,(array,title) in enumerate([(syn,'syn state'),(unstack,'unstacked from sequence neighbors')]):
            axes[row,col].imshow(array,aspect='auto',origin='lower',extent=[0,5,1,40],cmap='viridis',vmin=0,vmax=1,interpolation='nearest')
            axes[row,col].set_title(f'{label}: {title}');axes[row,col].set_xlabel('Window time (ns)')
    axes[0,0].set_ylabel('Residue');axes[1,0].set_ylabel('Residue');save_figure(fig,out,'03_base_state_heatmaps')

    fig,axes=plt.subplots(2,2,figsize=(10,9),layout='constrained')
    for col,(label,result) in enumerate(results.items()):
        im=axes[0,col].imshow(result['contact'],origin='lower',extent=[.5,40.5,.5,40.5],vmin=0,vmax=1,cmap='magma');axes[0,col].set_title(label+' nonlocal contacts')
        axes[1,col].imshow(result['stack'],origin='lower',extent=[.5,40.5,.5,40.5],vmin=0,vmax=1,cmap='magma');axes[1,col].set_title(label+' stacking')
        for ax in axes[:,col]:ax.set_xlabel('Residue');ax.set_ylabel('Residue')
    fig.colorbar(im,ax=axes.ravel().tolist(),label='Frame occupancy',shrink=.75,pad=.02);save_figure(fig,out,'04_contact_and_stack_maps')

    fig,axes=plt.subplots(2,2,figsize=(10,7))
    for ax,(key,title) in zip(axes.flat,[('rg_A','Rg (Å)'),('end_to_end_A',"5′–3′ distance (Å)"),('adjacent_stacks','Adjacent stacks'),('syn_fraction','syn fraction')]):
        for label,result in results.items():ax.hist(result['frame'][key],bins=25,density=True,histtype='step',lw=2,color=COLORS[label],label=label)
        ax.set_xlabel(title);ax.set_ylabel('Density');ax.grid(alpha=.2)
    axes[0,0].legend(frameon=False);save_figure(fig,out,'05_distributions')

    fig,axes=plt.subplots(3,1,figsize=(10,8),sharex=True)
    for ax,key,title in zip(axes,['K_within_3p5A','Cl_within_3p5A','waterO_within_3p5A'],['K⁺ within 3.5 Å','Cl⁻ within 3.5 Å','Water O within 3.5 Å']):
        for label,result in results.items():ax.plot(result['environment'].window_time_ns,result['environment'][key],'-o',ms=2,lw=1.2,color=COLORS[label],label=label)
        ax.set_ylabel(title);ax.grid(alpha=.2)
    axes[-1].set_xlabel('Time in analyzed window (ns)');axes[0].legend(frameon=False);save_figure(fig,out,'06_ion_and_hydration_contacts')

    fig,axes=plt.subplots(2,2,figsize=(10,9),layout='constrained')
    for row,(label,result) in enumerate(results.items()):
        rep=result['coords'][result['validation']['representative_frame_zero_based']-result['validation']['selected_frame_range_zero_based'][0]]
        c=rep[result['c1']];bc=result['base_centers'][result['validation']['representative_frame_zero_based']-result['validation']['selected_frame_range_zero_based'][0]]
        for col,(a,b,names) in enumerate([(0,2,('x','z')),(0,1,('x','y'))]):
            ax=axes[row,col];ax.plot(c[:,a],c[:,b],'-',color=COLORS[label],lw=1.5);points=ax.scatter(bc[:,a],bc[:,b],c=np.arange(1,41),cmap='viridis',s=18)
            ax.set_aspect('equal',adjustable='datalim');ax.set_xlabel(names[0]+' (Å)');ax.set_ylabel(names[1]+' (Å)');ax.set_title(label+' representative')
    fig.colorbar(points,ax=axes.ravel().tolist(),label='Residue',shrink=.7,pad=.02);save_figure(fig,out,'07_representative_structures')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--da-psf',type=Path,required=True);parser.add_argument('--da-dcd',type=Path,required=True)
    parser.add_argument('--dt-psf',type=Path,required=True);parser.add_argument('--dt-dcd',type=Path,required=True)
    parser.add_argument('--frames',type=int,default=500)
    parser.add_argument('--da-start-frame',type=int,default=0)
    parser.add_argument('--dt-start-frame',type=int,default=0)
    parser.add_argument('--out',type=Path,default=Path(__file__).resolve().parent.parent/'results')
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    setups={'dA40':source_setup('dA40',args.da_psf,args.da_dcd,args.frames,args.da_start_frame),
            'dT40':source_setup('dT40',args.dt_psf,args.dt_dcd,args.frames,args.dt_start_frame)}
    results={label:analyze_one(info,args.out) for label,info in setups.items()}
    pd.concat([result['frame'] for result in results.values()]).to_csv(args.out/'timeseries.csv',index=False,float_format='%.8g')
    pd.concat([result['local'] for result in results.values()]).to_csv(args.out/'per_residue_timeseries.csv',index=False,float_format='%.8g')
    pd.concat([result['summary'] for result in results.values()]).to_csv(args.out/'per_residue_summary.csv',index=False,float_format='%.8g')
    pd.concat([result['environment'] for result in results.values()]).to_csv(args.out/'environment_timeseries.csv',index=False,float_format='%.8g')
    summaries=[]
    for label,result in results.items():
        frame=result['frame'];env=result['environment']
        row={'sequence':label}
        for key in ['rg_A','end_to_end_A','c1_contour_A','extension_fraction','shape_anisotropy','adjacent_stacks','nonadjacent_stacks','unstacked_neighbor_fraction','unstacked_any_fraction','syn_fraction','outward_projection_mean','nonlocal_contact_pairs','base_base_hbonds','all_interresidue_hbonds']:
            row[key+'_mean']=frame[key].mean();row[key+'_sd']=frame[key].std(ddof=1)
        for key in ['K_within_3p5A','Cl_within_3p5A','waterO_within_3p5A','K_box_molar','Cl_box_molar']:
            row[key+'_mean']=env[key].mean();row[key+'_sd']=env[key].std(ddof=1)
        summaries.append(row)
        for name in ['contact','stack','hbond']:
            pd.DataFrame(result[name],index=range(1,41),columns=range(1,41)).to_csv(args.out/f'{label}_{name}_occupancy.csv')
    summary=pd.DataFrame(summaries);summary.to_csv(args.out/'summary.csv',index=False,float_format='%.8g')
    blocks=pd.concat([r['frame'] for r in results.values()]);blocks['block_0p5ns']=np.minimum(9,np.floor((blocks.window_time_ns-1e-9)/.5).astype(int))
    numeric=[c for c in blocks.columns if c not in {'sequence','window_time_ns','source_time_ns','block_0p5ns'}]
    blocks.groupby(['sequence','block_0p5ns'])[numeric].mean().reset_index().to_csv(args.out/'block_means_0p5ns.csv',index=False,float_format='%.8g')
    validation={label:result['validation'] for label,result in results.items()}
    validation['method']='PSF-bond reconstruction; matched 0.01–5.00 ns production windows; frame statistics are descriptive'
    validation['input_sha256']={
        'dA40_system.psf':sha256(args.da_psf),
        'dA40_production.dcd':sha256(args.da_dcd),
        'dT40_system.psf':sha256(args.dt_psf),
        'dT40_production.dcd':sha256(args.dt_dcd),
    }
    (args.out/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    plot_all(results,args.out)
    print(summary.to_string(index=False),flush=True)


if __name__=='__main__':main()
