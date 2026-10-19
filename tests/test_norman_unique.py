"""Duplicate-symbol preparation preserves integer counts, cells and source identity."""
import json
from pathlib import Path
import sys
import tempfile

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from prep_norman_unique import prepare


def main():
    with tempfile.TemporaryDirectory() as directory:
        inp, out = Path(directory)/"in.h5ad", Path(directory)/"out.h5ad"
        counts = np.array([[1,2,3],[4,5,6],[7,8,9]], dtype=np.float32)
        a = ad.AnnData(sparse.csr_matrix(counts),
            obs=pd.DataFrame({"perturbation":["control","A","B"],"gemgroup":[1,1,2]},index=["c0","c1","c2"]),
            var=pd.DataFrame({"ensembl_id":["e1","e2","e3"]},index=["A","A","B"]))
        a.write_h5ad(inp)
        original = inp.read_bytes()
        note = prepare(inp,out)
        b = ad.read_h5ad(out)
        assert b.var_names.tolist()==["A","B"]
        np.testing.assert_array_equal(b.X.toarray(),[[3,3],[9,6],[15,9]])
        np.testing.assert_array_equal(b.X.sum(1),a.X.sum(1))
        pd.testing.assert_frame_equal(b.obs,a.obs)
        assert b.var['source_ensembl_ids'].tolist()==["e1;e2","e3"]
        assert note['duplicate_columns_aggregated']==1 and note['per_cell_library_totals_exactly_preserved']
        assert inp.read_bytes()==original
        assert json.loads(out.with_suffix('.manifest.json').read_text())==note
        before=out.read_bytes()
        try: prepare(inp,out)
        except FileExistsError: pass
        else: raise AssertionError('existing output was overwritten')
        assert out.read_bytes()==before
        a.X=sparse.csr_matrix(counts/10)
        a.write_h5ad(inp)
        try: prepare(inp,Path(directory)/'invalid.h5ad')
        except ValueError: pass
        else: raise AssertionError('nonraw counts were accepted')
    print('ALL NORMAN-UNIQUE TESTS PASSED')


if __name__=='__main__': main()
