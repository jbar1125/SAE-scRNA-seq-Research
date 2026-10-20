"""Both component signs and discovery-only selection need meaningful checks."""
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from norman_direct_pca import analyze


def main():
    rng=np.random.default_rng(197)
    n=400;labels=np.array(['control']*200+['P0']*100+['P1']*100)
    z=rng.normal(size=(n,2));z[labels=='P0',0]-=3;z[labels=='P1',1]+=3
    x=rng.normal(size=(n,100));x[:,:5]=-z[:,0,None]+rng.normal(0,.02,(n,5))
    x[:,5:10]=z[:,1,None]+rng.normal(0,.02,(n,5))
    genes=[f'g{i}' for i in range(100)]
    panels={'negative_axis':genes[:5],'positive_axis':genes[5:10],'absent':['missing1','missing2','missing3']}
    strata=np.array(['batch1','batch2']*(n//2))
    a=analyze(z,x,genes,labels,panels,strata,top_k=10,n_shuffle=99,min_active=20)
    assert a['family_size']==6
    assert a['programs']['negative_axis']['modulators']==['P0']
    assert a['programs']['negative_axis']['orientation']=='negative'
    assert a['programs']['positive_axis']['modulators']==['P1']
    assert a['programs']['absent']['low_power']
    assert a['identity_heterogeneity']['programs']['absent']['empirical_p']==1
    b=analyze(-z,x,genes,labels,panels,strata,top_k=10,n_shuffle=99,min_active=20)
    for name in panels:
        assert a['programs'][name]['modulators']==b['programs'][name]['modulators']
        assert a['programs'][name]['enrich_p_bonf']==b['programs'][name]['enrich_p_bonf']
    changed=x.copy();changed[a['control_split']['inference_rows']]=rng.normal(0,100,(100,100))
    c=analyze(z,changed,genes,labels,panels,strata,top_k=10,n_shuffle=99,min_active=20)
    for name in panels:
        assert a['programs'][name]['feature']==c['programs'][name]['feature']
        assert a['programs'][name]['enrich_p_bonf']==c['programs'][name]['enrich_p_bonf']
    print('ALL DIRECT-PCA BASELINE TESTS PASSED')


if __name__=='__main__':main()
