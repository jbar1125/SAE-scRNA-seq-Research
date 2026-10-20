"""Partial or heterogeneous runs must never become a completed three-seed result."""
import copy
from pathlib import Path
import sys
import json
import tempfile
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from archive_norman_cpu_results import PANELS,INPUT_SHA,validate_record,joint_summary
import archive_norman_cpu_results as archive


def fixture(seed):
    return {'sae':{'epochs':200,'epochs_completed':200,'latent_dim':2048,'k':32,
        'batch':1024,'trained_on':'all','n_train_cells':111445,'d':4081},
        'provenance':{'arguments':{'seed':seed,'strata_col':'gemgroup'},
            'input_sha256':INPUT_SHA,'source_sha256':{'fixed_source':'a'}},
        'programs':{name:{'low_power':i!=0,'n_modulators':2 if i==0 else None,
            'modulators':['shared','only'+str(seed)] if i==0 else []} for i,name in enumerate(PANELS)},
        'identity_heterogeneity':{'n_planned_programs':9,'n_permutations_requested':999,
            'n_permutations':999,'stratified':True,'programs':{
                name:{'matched':i==0,'empirical_p':.001 if i==0 else 1.} for i,name in enumerate(PANELS)}}}


def refuses(call):
    try:call()
    except ValueError:return
    raise AssertionError('invalid completion or aggregation was accepted')


def main():
    r=fixture(0);validate_record(r,0)
    for key,value in [('epochs_completed',20),('latent_dim',512),('trained_on','control'),('d',4000)]:
        bad=copy.deepcopy(r);bad['sae'][key]=value
        refuses(lambda:validate_record(bad,0))
    bad=copy.deepcopy(r);bad['identity_heterogeneity']['programs'].pop(PANELS[-1])
    refuses(lambda:validate_record(bad,0))
    bad=copy.deepcopy(r);bad['identity_heterogeneity']['n_permutations']=0
    refuses(lambda:validate_record(bad,0))
    bad=copy.deepcopy(r);bad['identity_heterogeneity']['programs'][PANELS[0]]['empirical_p']=0.
    refuses(lambda:validate_record(bad,0))
    rows=[fixture(s) for s in range(3)]
    result=joint_summary(rows)
    assert result['family_size']==27 and result['complete_seeds']==[0,1,2]
    p=result['programs'][PANELS[0]]
    assert p['joint_identity_rejections']==3 and p['matched_seeds']==3
    assert all(np.isclose(v,.027) for v in p['joint_holm_p_by_seed'].values())
    assert p['descriptive_candidates_in_at_least_two_seeds']==['shared']
    assert result['programs'][PANELS[1]]['matched_seeds']==0
    assert result['programs'][PANELS[1]]['joint_identity_rejections']==0
    refuses(lambda:joint_summary(rows[:2]))
    refuses(lambda:joint_summary([fixture(0),fixture(0),fixture(2)]))
    rows[2]['provenance']['source_sha256']['fixed_source']='changed'
    refuses(lambda:joint_summary(rows))
    with tempfile.TemporaryDirectory() as tmp:
        root=Path(tmp);run=root/'run';dest=root/'archive';run.mkdir()
        source=run/'seed0.json';source.write_text(json.dumps(fixture(0)))
        import os
        os.utime(source,(0,0))
        (run/'seed0.log').write_text('completed training\n')
        with patch.object(archive,'ROOT',root),patch.object(archive,'RUN',run),patch.object(archive,'ARCHIVE',dest),patch.object(archive,'verify_large_artifacts') as verify,patch.object(archive,'publish') as publish:
            assert not archive.collect()
            verify.assert_called_once()
            original=(dest/'seed0.json').read_bytes()
            # Recover a copy completed before manifest creation/publication.
            (dest/'manifest.json').unlink()
            assert not archive.collect(publish_git=True)
            assert (dest/'seed0.json').read_bytes()==original
            assert (dest/'manifest.json').exists()
            assert dest/'seed0.json' in publish.call_args.args[0]
            assert not archive.collect(publish_git=True)
            assert verify.call_count==1
        existing=dest/'seed0.json'
        try:archive.copy_exclusive(source,existing)
        except FileExistsError:pass
        else:raise AssertionError('archive copy overwrote an existing file')
        assert existing.read_bytes()==original
        assert sorted(p.name for p in dest.iterdir())==['manifest.json','seed0.json','seed0.log']
    print('ALL NORMAN-RESULT ARCHIVE TESTS PASSED')


if __name__=='__main__':main()
