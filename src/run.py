import numpy as np,pandas as pd,json
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score,f1_score
from sklearn.preprocessing import OneHotEncoder
P='/tmp/poore/'; X=np.load(P+'X.npy'); idx=pd.read_csv(P+'idx.csv').iloc[:,0].tolist(); cols=pd.read_csv(P+'cols.csv').iloc[:,0].tolist()
m=pd.read_csv(P+'meta.csv',index_col=0).loc[idx].reset_index(drop=True)
early=m.pathologic_stage_label.str.fullmatch(r'Stage (I|IA|IB|IS|II|IIA|IIB|IIC)|I or II NOS'); late=m.pathologic_stage_label.str.fullmatch(r'Stage (III|IIIA|IIIB|IIIC|IV|IVA|IVB|IVC)')
T=m.sample_type; PT=T=='Primary Tumor'; STN=T=='Solid Tissue Normal'; BN=T=='Blood Derived Normal'
R={}
R['d1']=dict(n_samples=len(m),n_genera=X.shape[1],early_primary=int((PT&early).sum()),late_primary=int((PT&late).sum()),solid_normal=int(STN.sum()),blood_normal=int(BN.sum()),
  solid_normal_by_type=m[STN].disease_type.value_counts().to_dict(),early_by_type_with_normals={k:int(v) for k,v in m[PT&early&m.disease_type.isin(m[STN].disease_type.unique())].disease_type.value_counts().items()})
enc=OneHotEncoder(handle_unknown='ignore').fit(m[['disease_type','data_submitting_center_label','experimental_strategy']])
MZ=enc.transform(m[['disease_type','data_submitting_center_label','experimental_strategy']])
enc2=OneHotEncoder(handle_unknown='ignore').fit(m[['data_submitting_center_label','experimental_strategy']]); MZ2=enc2.transform(m[['data_submitting_center_label','experimental_strategy']])
mu=X.mean(0); sd=X.std(0)+1e-6; Xs=((X-mu)/sd).astype(np.float32)
def lr(Xa,y,C=.1): return LogisticRegression(C=C,max_iter=500).fit(Xa,y)
def cv(sel,y,feat='mic',k=5,shuffle=False,coef=False):
    ii=np.where(sel)[0]; yy=y[ii]; g=m.case_uuid.values[ii]; k=min(k,len(set(g))); pred=np.zeros(len(ii)); tops=[]
    rng=np.random.default_rng(0)
    for tr,te in GroupKFold(k).split(ii,yy,g):
        ytr=rng.permutation(yy[tr]) if shuffle else yy[tr]
        if len(set(ytr))<2: pred[te]=.5; continue
        A=Xs[ii] if feat=='mic' else (MZ[ii] if feat=='meta' else MZ2[ii])
        mdl=lr(A[tr],ytr); pred[te]=mdl.decision_function(A[te])
        if coef: tops.append(set(np.argsort(-np.abs(mdl.coef_[0]))[:20]))
    a=roc_auc_score(yy,pred) if len(set(yy))>1 else np.nan
    return (a,tops) if coef else a
sel=(PT&early&True)|STN; sel=sel.values; y=(PT&early).values.astype(int)
sel2=sel&(m.disease_type.isin(m[STN].disease_type.unique())).values
a_mic,tops=cv(sel2,y,'mic',coef=True); a_meta=cv(sel2,y,'meta'); a_meta2=cv(sel2,y,'meta2')
R['d2']=dict(n=int(sel2.sum()),n_tumor=int(y[sel2].sum()),n_normal=int((1-y[sel2]).sum()),microbial_auroc=round(a_mic,3),metadata_disease_center_strategy_auroc=round(a_meta,3),metadata_center_strategy_only_auroc=round(a_meta2,3))
R['d2_gate']=bool(a_mic>=.80 and a_mic>=a_meta+.05); print(R['d2'],flush=True)
J=[len(a&b)/len(a|b) for i,a in enumerate(tops) for b in tops[i+1:]]
cnt=pd.Series([cols[i] for t in tops for i in t]).value_counts()
R['d9']=dict(mean_pairwise_jaccard_top20=round(float(np.mean(J)),3),genera_in_all_5_folds=cnt[cnt==5].index.tolist()[:20])
R['d7_shuffle_auroc']=round(float(np.mean([cv(sel2,y,'mic',shuffle=True)])),3)
# D3 and D4
types=[t for t in m[STN].disease_type.unique() if (STN&(m.disease_type==t)).sum()>=15 and (PT&early&(m.disease_type==t)).sum()>=15]
d3={};d4={}
for t in types:
    s=sel&(m.disease_type==t).values
    d3[t]=dict(n_tumor=int(y[s].sum()),n_normal=int((1-y[s]).sum()),microbial=round(cv(s,y,'mic'),3),meta_center_strategy=round(cv(s,y,'meta2'),3))
    tr=sel2&(m.disease_type!=t).values; te=s; mdl=lr(Xs[tr],y[tr]); d4[t]=round(float(roc_auc_score(y[te],mdl.decision_function(Xs[te]))),3)
    print(t,d3[t],d4[t],flush=True)
R['d3']=d3; R['d4']=dict(per_type=d4,mean=round(float(np.mean(list(d4.values()))),3)); R['d4_gate_mean_ge_0.70']=bool(R['d4']['mean']>=.70)
# D5
nrm=STN.values&sel2; par=(pd.util.hash_array(m.case_uuid.values)%2==0)
trn=((PT&late).values&m.disease_type.isin(m[STN].disease_type.unique()).values)|(nrm&par); ten=((PT&early).values&m.disease_type.isin(m[STN].disease_type.unique()).values)|(nrm&~par)
ytr=(PT.values[trn]).astype(int); yte=PT.values[ten].astype(int); mdl=lr(Xs[trn],ytr); mm=lr(MZ[trn],ytr)
R['d5']=dict(n_train=int(trn.sum()),n_test=int(ten.sum()),microbial=round(float(roc_auc_score(yte,mdl.decision_function(Xs[ten]))),3),metadata=round(float(roc_auc_score(yte,mm.decision_function(MZ[ten]))),3))
# D6
s6=(PT&early).values; yc=m.disease_type.values; ii=np.where(s6)[0]; pr=np.empty(len(ii),dtype=object); prm=np.empty(len(ii),dtype=object)
for tr,te in GroupKFold(5).split(ii,yc[ii],m.case_uuid.values[ii]):
    pr[te]=LogisticRegression(C=.1,max_iter=300).fit(Xs[ii][tr],yc[ii][tr]).predict(Xs[ii][te]); prm[te]=LogisticRegression(max_iter=300).fit(MZ2[ii][tr],yc[ii][tr]).predict(MZ2[ii][te])
R['d6']=dict(n=len(ii),n_classes=len(set(yc[ii])),macroF1_microbial=round(float(f1_score(yc[ii],pr,average='macro')),3),macroF1_center_strategy_only=round(float(f1_score(yc[ii],prm,average='macro')),3),majority_acc=round(float(pd.Series(yc[ii]).value_counts(normalize=True).iloc[0]),3))
# D10
s10=((PT&early)|BN).values; R['d10']=dict(n_tumor=int((PT&early).sum()),n_blood_normal=int(BN.sum()),microbial=round(cv(s10,y,'mic'),3),meta=round(cv(s10,y,'meta'),3),note='confounded: different tissue and sample type')
json.dump(R,open('results/results.json','w'),indent=1); print(json.dumps(R,indent=1)[:6000])
