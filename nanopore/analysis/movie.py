"""CPU ray-traced VMD frames -> H.264 MP4, annotated with measured quantities."""
import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from PIL import Image,ImageDraw,ImageFont
import imageio_ffmpeg

ROOT=Path(__file__).resolve().parent
FONT='/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'
BOLD='/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'


def compose(picture,row,summary):
    width,height=1280,960
    canvas=Image.new('RGB',(width,height),'white')
    picture=picture.convert('RGB');picture.thumbnail((950,830),Image.Resampling.LANCZOS)
    canvas.paste(picture,((950-picture.width)//2,90+(830-picture.height)//2))
    draw=ImageDraw.Draw(canvas)
    def text(x,y,s,size=22,color='#24384a',bold=False):
        draw.text((x,y),s,font=ImageFont.truetype(BOLD if bold else FONT,size),fill=color)
    draw.rectangle((0,0,width,80),fill='#153348')
    text(28,17,f"poly(dA)40 / {summary['pore'].upper()}   |   nominal 3 M KCl",28,'white',True)
    status='5 ns production segment' if summary['status']=='complete' else 'IN-PROGRESS PREVIEW'
    text(28,51,status,17,'#cae2ee')
    draw.line((959,110,959,875),fill='#dfe6ea',width=2)
    text(985,116,'SIMULATION TIME',16,'#647a89',True)
    text(985,146,f"{row['time_ns']:.3f} ns",32,bold=True)
    text(985,189,f"of {summary['production_ns']:.3f} ns shown",17)
    text(985,257,'DNA conformation',19,bold=True)
    text(985,294,f"Rg: {row['rg_A']:.2f} A",21)
    text(985,329,f"End-to-end: {row['end_to_end_A']:.1f} A",18)
    text(985,390,'Protein contacts',19,bold=True)
    text(985,428,f"{int(row['contacting_DNA_residues_3p5A'])} / 40 nucleotides",21)
    text(985,465,'Heavy-atom cutoff: 3.5 A',16)
    text(985,523,'Axial advance of 5\' end',17,bold=True)
    advance=summary['initial']['lead_z_A']-row['lead_z_A']
    text(985,557,f"{advance:+.2f} A toward -z",19)
    for y,color,label in [(633,'#f28f1f','DNA'),(670,'#518cb5','Pore protein'),(707,'#d62e29',"5' nucleotide"),(744,'#30a663',"3' nucleotide"),(781,'#a1a6ad','Membrane P atoms')]:
        draw.ellipse((987,y+3,1000,y+16),fill=color);text(1011,y,label,17)
    text(28,903,'Actual saved coordinates | 25 ps/frame | water and ions hidden | front protein chains transparent',16,'#566b7b')
    text(28,930,'Pore-aligned view. Axial motion and contacts do not by themselves establish complete passage.',15,'#566b7b')
    end=summary['production_ns']
    draw.rectangle((0,82,int(width*row['time_ns']/end),87),fill='#ee962c')
    return canvas


def main():
    ap=argparse.ArgumentParser();ap.add_argument('input',type=Path)
    ap.add_argument('--test-frame',type=int);ap.add_argument('--encode-only',action='store_true')
    a=ap.parse_args();folder=a.input.resolve()
    config=json.loads((folder/'render.json').read_text());summary=json.loads((folder/'summary.json').read_text())
    table=pd.read_csv(folder/'timeseries.csv');frames=folder/'raytrace';frames.mkdir(exist_ok=True)
    first=0 if a.test_frame is None else a.test_frame;last=len(table)-1 if a.test_frame is None else a.test_frame
    if not a.encode_only:
        env=dict(os.environ,VMDFORCECPUCOUNT='2',VMDCUDA='0',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
        span=max(config['z_bounds'][1]-config['z_bounds'][0],150.)
        command=[os.environ.get('VMD_BIN','/opt/bin/vmd'),'-dispdev','text','-e',str(ROOT/'render.tcl'),'-args',str(folder),str(frames),
                 *map(str,config['center']),str(span),str(first),str(last)]
        log=folder/('render_test.log' if a.test_frame is not None else 'render.log')
        with log.open('w') as f:
            # Keep stdin open until the Tcl script explicitly quits. The old
            # headless VMD build can tear down its graphics on inherited EOF.
            proc=subprocess.Popen(command,env=env,stdin=subprocess.PIPE,stdout=f,stderr=subprocess.STDOUT)
            rc=proc.wait()
            proc.stdin.close()
        data=log.read_text(errors='replace')
        if rc!=0 or 'RENDER_COMPLETE' not in data or 'RENDER_FAILED' in data:
            raise RuntimeError('VMD rendering failed; inspect '+str(log))
        checks=pd.read_csv(folder/f'render_checks_{first}_{last}.csv')
        for _,check in checks.iterrows():
            source=table.iloc[int(check['frame'])]
            for key in ['rg_A','end_to_end_A']:
                if not np.isclose(source[key],check[key],atol=1e-3,rtol=1e-5):
                    raise ValueError('Independent VMD measurement disagrees: '+key)
            if int(source['contact_pairs_3p5A'])!=int(check['contact_pairs_3p5A']):
                raise ValueError('Independent VMD contact count disagrees')
    if a.test_frame is not None:
        compose(Image.open(frames/f'frame_{first:04d}.tga'),table.iloc[first],summary).save(folder/'render_test.png')
        print(folder/'render_test.png');return
    movie=folder/'trajectory.mp4'
    writer=imageio_ffmpeg.write_frames(str(movie),(1280,960),fps=10,codec='libx264',pix_fmt_out='yuv420p',
                                      quality=8,macro_block_size=16,output_params=['-movflags','+faststart','-threads','2'])
    writer.send(None)
    distinct=set()
    try:
        for i,row in table.iterrows():
            picture=Image.open(frames/f'frame_{i:04d}.tga').convert('RGB')
            pixels=np.asarray(picture)
            if np.mean(pixels.min(axis=2)<220)<.005: raise ValueError('Blank molecular rendering')
            distinct.add(hashlib.sha256(pixels.tobytes()).hexdigest())
            canvas=compose(picture,row,summary)
            if i==len(table)-1: canvas.save(folder/'poster.png')
            writer.send(np.asarray(canvas))
    finally: writer.close()
    if len(distinct)<len(table)*.9: raise ValueError('Frozen or repeated molecular animation frames')
    # Decode the complete MP4 back to frames; metadata/file existence alone is
    # not enough to verify that it contains every requested trajectory frame.
    reader=imageio_ffmpeg.read_frames(str(movie),pix_fmt='rgb24')
    meta=next(reader);count=sum(1 for _ in reader)
    if count!=len(table): raise ValueError(f'Encoded {count} frames, expected {len(table)}')
    result=dict(codec='H.264',container='MP4',dimensions=[1280,960],fps=10,frames=count,
                duration_seconds=count/10,physical_time_ns=summary['production_ns'],source_frame_spacing_ps=25,
                interpolation=False,decoded_all_frames=True,ffmpeg_version=imageio_ffmpeg.get_ffmpeg_version(),
                distinct_molecular_frames=len(distinct),water_and_ions_hidden=True,protein_front_chains_transparent=True)
    (folder/'movie_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__': main()
