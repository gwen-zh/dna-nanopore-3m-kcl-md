"""Small PSF/NAMD DCD reader and display-only DCD writer (orthorhombic)."""
import collections
import struct
from pathlib import Path
import numpy as np
import pandas as pd


def psf(path):
    with Path(path).open() as f:
        for line in f:
            if '!NATOM' in line: n=int(line.split()[0]); break
        rows=[next(f).split() for _ in range(n)]
        for line in f:
            if '!NBOND' in line: nb=int(line.split()[0]); break
        ids=[]
        while len(ids)<2*nb: ids.extend(map(int,next(f).split()))
    atoms=pd.DataFrame(rows,columns=['id','segid','resid','resname','name','type','charge','mass','unused'])
    for name in ['id','resid']: atoms[name]=atoms[name].astype(int)
    for name in ['charge','mass']: atoms[name]=atoms[name].astype(float)
    return atoms,np.array(ids).reshape(-1,2)-1


def edges(n,bonds):
    adjacent=[[] for _ in range(n)]
    for i,j in bonds: adjacent[i].append(j); adjacent[j].append(i)
    seen={0}; q=collections.deque([0]); result=[]
    while q:
        i=q.popleft()
        for j in adjacent[i]:
            if j not in seen: seen.add(j);q.append(j);result.append((i,j))
    if len(seen)!=n: raise ValueError('Disconnected molecular fragment')
    return result


def whole(x,box,tree):
    joined=x.copy()
    for i,j in tree:
        d=x[j]-x[i];joined[j]=joined[i]+d-box*np.rint(d/box)
    return joined


def record(f):
    b=f.read(4)
    if len(b)!=4: raise EOFError('Incomplete DCD record')
    n=struct.unpack('<i',b)[0]
    if not 0<=n<100000000: raise ValueError('Bad DCD record size')
    data=f.read(n)
    if len(data)!=n or f.read(4)!=b: raise EOFError('Incomplete DCD payload/footer')
    return data


class DCD:
    def __init__(self,path,n,complete):
        self.path=Path(path);self.f=self.path.open('rb');self.header=record(self.f)
        if len(self.header)!=84 or self.header[:4]!=b'CORD': raise ValueError('Unsupported DCD')
        c=np.frombuffer(self.header[4:],dtype='<i4')
        self.first,self.stride=int(c[1]),int(c[2])
        if c[8]!=0 or c[10]!=1: raise ValueError('Fixed atoms or nonperiodic source unsupported')
        record(self.f)
        if struct.unpack('<i',record(self.f))[0]!=n: raise ValueError('DCD/PSF atom count mismatch')
        self.n=n;self.offset=self.f.tell();self.framebytes=56+3*(8+4*n)
        payload=self.path.stat().st_size-self.offset
        self.count,self.trailing=divmod(payload,self.framebytes)
        if complete and (self.trailing or self.count!=int(c[0])):
            self.f.close()
            raise ValueError('Incomplete closed DCD')
        self.byte_limit=self.offset+self.count*self.framebytes

    def __iter__(self):
        for index in range(self.count):
            self.f.seek(self.offset+index*self.framebytes)
            box=np.frombuffer(record(self.f),dtype='<f8')
            if len(box)!=6 or not np.all((abs(box[[1,3,4]])<1e-6)|(abs(box[[1,3,4]]-90)<1e-6)):
                raise ValueError('Nonorthorhombic source box')
            xyz=np.column_stack([np.frombuffer(record(self.f),dtype='<f4') for _ in range(3)]).astype(float)
            if xyz.shape!=(self.n,3) or not np.isfinite(xyz).all(): raise ValueError('Invalid source frame')
            yield self.first+index*self.stride,xyz,box[[0,2,5]]


def binary(path,n):
    with Path(path).open('rb') as f:
        if struct.unpack('<i',f.read(4))[0]!=n: raise ValueError('Binary atom count mismatch')
        xyz=np.fromfile(f,dtype='<f8')
    if xyz.size!=n*3 or not np.isfinite(xyz).all(): raise ValueError('Invalid binary coordinates')
    return xyz.reshape(n,3)


def write_record(f,data):
    f.write(struct.pack('<i',len(data)));f.write(data);f.write(struct.pack('<i',len(data)))


def write_dcd(path,xyz,header):
    # Display coordinates are rigidly aligned to the pore. Do NOT attach the
    # unrotated periodic cell: these are visualization-only, not restart data.
    control=np.frombuffer(header[4:],dtype='<i4').copy()
    control[0]=len(xyz);control[1]=0;control[2]=25000;control[3]=(len(xyz)-1)*25000
    control[10]=0
    with Path(path).open('xb') as f:
        write_record(f,b'CORD'+control.tobytes())
        write_record(f,struct.pack('<i',1)+b'Pore-aligned display trajectory; not an MD restart'.ljust(80))
        write_record(f,struct.pack('<i',xyz.shape[1]))
        for frame in xyz:
            for axis in range(3): write_record(f,np.asarray(frame[:,axis],dtype='<f4').tobytes())


def write_psf(path,atoms,bonds):
    # Preserve the selected atom identities and bonds, without pretending that
    # omitted solvent/force-field interaction terms form an MD-ready system.
    with Path(path).open('x') as f:
        f.write('PSF\n\n       1 !NTITLE\n REMARKS Visualization-only subset, not MD topology\n\n')
        f.write(f'{len(atoms):8d} !NATOM\n')
        for i,(_,a) in enumerate(atoms.iterrows(),1):
            f.write(f'{i:8d} {a.segid:<4} {a.resid:<4d} {a.resname:<4} {a["name"]:<4} {a.type:<6} {a.charge:10.6f} {a.mass:10.4f}           0\n')
        f.write(f'\n{len(bonds):8d} !NBOND: bonds\n')
        values=(bonds+1).ravel()
        for start in range(0,len(values),8): f.write(''.join(f'{v:8d}' for v in values[start:start+8])+'\n')
        for title in ['NTHETA: angles','NPHI: dihedrals','NIMPHI: impropers','NDON: donors','NACC: acceptors']:
            f.write(f'\n       0 !{title}\n')
        f.write('\n       0 !NNB\n\n')
        for start in range(0,len(atoms),8): f.write(''.join(f'{0:8d}' for _ in range(min(8,len(atoms)-start)))+'\n')
        f.write('\n       1       0 !NGRP NST2\n       0       0       0\n')


def write_pdb(path,atoms,xyz):
    with Path(path).open('x') as f:
        f.write('REMARK Visualization-only protein-aligned subset; solvent omitted\n')
        for i,((_,a),p) in enumerate(zip(atoms.iterrows(),xyz),1):
            name=a['name']; element='H' if a.mass<2 else name[0]
            chain='D' if a.segid=='AN1' else a.segid[0]
            f.write(f'ATOM  {i:5d} {name:^4} {a.resname:<4}{chain}{a.resid:4d}    '
                    f'{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}{1.:6.2f}{0.:6.2f}      {a.segid:<4}{element:>2}\n')
        f.write('END\n')
