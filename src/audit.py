import numpy as np,pandas as pd,json
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import OneHotEncoder
P='/tmp/poore/'
v=pd.read_csv(P+'voom.csv',index_col=0); X=v.values.astype(np.float32); idx=v.index.tolist(); del v
m=pd.read_csv(P+'meta.csv',index_col=0).loc[idx].reset_index(drop=True)
early=m.pathologic_stage_label.str.fullmatch(r'Stage (I|IA|IB|IS|II|IIA|IIB|IIC)|I or II NOS'); T=m.sample_type; PT=T=='Primary Tumor'; STN=T=='Solid Tissue Normal'
sel=(((PT&early)|STN)&m.disease_type.isin(m[STN].disease_type.unique())).values; y=(PT&early).values.astype(int)
ii=np.where(sel)[0]; print('n',len(ii),y[ii].sum(),flush=True)
Xs=((X-X.mean(0))/(X.std(0)+1e-6)).astype(np.float32)[ii]; yy=y[ii]; mm=m.iloc[ii].reset_index(drop=True); g=mm.case_uuid.values
lr=lambda A,t: LogisticRegression(C=.1,max_iter=500).fit(A,t)
R={}
# D2 reproduction + Q2 predictions
pred=np.zeros(len(ii))
for tr,te in GroupKFold(5).split(ii,yy,g): pred[te]=lr(Xs[tr],yy[tr]).decision_function(Xs[te])
R['D2_repro_auroc']=float(roc_auc_score(yy,pred)); print('D2',R['D2_repro_auroc'],flush=True)
# Q2
st=(mm.disease_type+'|'+mm.data_submitting_center_label+'|'+mm.experimental_strategy).values; rows=[]
for s in sorted(set(st)):
    k=st==s
    if yy[k].sum()>=15 and (1-yy[k]).sum()>=15: rows.append((s,int(k.sum()),int(yy[k].sum()),float(roc_auc_score(yy[k],pred[k]))))
w=np.array([r[1] for r in rows]); a=np.array([r[3] for r in rows])
R['Q2']=dict(strata=rows,n_strata=len(rows),weighted_mean_auroc=float((w*a).sum()/w.sum()) if rows else None,min=float(a.min()) if rows else None); print('Q2',R['Q2']['n_strata'],R['Q2']['weighted_mean_auroc'],flush=True)
# Q1
cen=mm.data_submitting_center_label.values; q1=[]
for c in sorted(set(cen)):
    te=cen==c
    if yy[te].sum()>=20 and (1-yy[te]).sum()>=20:
        mdl=lr(Xs[~te],yy[~te]); q1.append((c,int(te.sum()),float(roc_auc_score(yy[te],mdl.decision_function(Xs[te])))))
R['Q1']=dict(centers=q1,mean_auroc=float(np.mean([r[2] for r in q1])) if q1 else None); print('Q1',R['Q1'],flush=True)
# Q3
cols=['disease_type','data_submitting_center_label','experimental_strategy','tissue_source_site_label','portion_ffpe','analyte_type_label','platform']
M=OneHotEncoder(handle_unknown='ignore').fit_transform(mm[cols].astype(str)); pm=np.zeros(len(ii))
for tr,te in GroupKFold(5).split(ii,yy,g): pm[te]=lr(M[tr],yy[tr]).decision_function(M[te])
R['Q3']=dict(extended_meta_auroc=float(roc_auc_score(yy,pm)),microbial_auroc=R['D2_repro_auroc']); R['Q3']['gate_met']=bool(R['D2_repro_auroc']>=.8 and R['D2_repro_auroc']>=R['Q3']['extended_meta_auroc']+.05); print('Q3',R['Q3'],flush=True)
json.dump(R,open('results/audit.json','w'),indent=1); print('done')
