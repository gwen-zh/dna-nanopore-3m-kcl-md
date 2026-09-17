import tempfile
import unittest
import struct
from pathlib import Path
import numpy as np
import pandas as pd
from trajectory import edges,whole,DCD,write_record,write_dcd,write_pdb
from analyze import rotation,sha


class GeometryAndIO(unittest.TestCase):
    def test_bonded_join(self):
        x=np.array([[9.5,0,0],[.5,0,0],[1.5,0,0]])
        y=whole(x,np.array([10.,10.,10.]),edges(3,np.array([[0,1],[1,2]])))
        np.testing.assert_allclose(y[:,0],[9.5,10.5,11.5])

    def test_disconnected_denied(self):
        with self.assertRaises(ValueError): edges(3,np.array([[0,1]]))

    def test_proper_rotation(self):
        x=np.random.default_rng(3).normal(size=(20,3))
        r=np.array([[0,-1,0],[1,0,0],[0,0,1]])
        target=x@r+[10,20,30]
        c,fit,ref=rotation(x,target)
        self.assertAlmostEqual(np.linalg.det(fit),1)
        np.testing.assert_allclose((x-c)@fit+ref,target,atol=1e-12)

    def source(self,path,truncated=False):
        c=np.zeros(20,dtype='<i4');c[0]=2;c[1]=25000;c[2]=25000;c[10]=1;c[19]=24
        with path.open('wb') as f:
            write_record(f,b'CORD'+c.tobytes());write_record(f,struct.pack('<i',1)+b'TEST'.ljust(80));write_record(f,struct.pack('<i',2))
            for k in range(2):
                write_record(f,np.array([10,0,10,0,0,10],dtype='<f8').tobytes())
                for axis in range(3): write_record(f,np.array([k+axis,k+axis+1],dtype='<f4').tobytes())
            if truncated: f.write(b'partial-tail')

    def test_growing_dcd_time_axis(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'source.dcd';self.source(p,True)
            reader=DCD(p,2,False)
            self.assertEqual([step for step,_,_ in reader],[25000,50000]);reader.f.close()
            with self.assertRaises(ValueError): DCD(p,2,True)

    def test_display_dcd_has_no_false_cell(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'source.dcd';self.source(p)
            r=DCD(p,2,True);xyz=np.zeros((3,2,3));out=Path(d)/'view.dcd'
            write_dcd(out,xyz,r.header);r.f.close()
            with out.open('rb') as f:
                f.read(8);c=np.frombuffer(f.read(80),dtype='<i4')
            self.assertEqual((c[0],c[1],c[10]),(3,0,0))

    def test_snapshot_hash(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x';p.write_bytes(b'headPAYLOADtail')
            before=sha(p,4,7);p.write_bytes(b'HEADPAYLOADlonger tail')
            self.assertEqual(before,sha(p,4,7))

    def test_pdb_columns(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x.pdb'
            a=pd.DataFrame([dict(segid='AN1',resid=40,resname='ADE',name="C1'",mass=12.)])
            write_pdb(p,a,np.array([[1.,2.,3.]]))
            line=p.read_text().splitlines()[1]
            self.assertEqual(line[12:16].strip(),"C1'");self.assertEqual(int(line[22:26]),40)
            self.assertEqual(line[72:76].strip(),'AN1')
            self.assertEqual([float(line[i:i+8]) for i in [30,38,46]],[1,2,3])


if __name__=='__main__': unittest.main()
