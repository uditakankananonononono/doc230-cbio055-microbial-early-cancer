"""Recovered preprocessing (verbatim logic from the original run): builds X.npy, idx.csv, cols.csv from voom.csv for src/run.py.
Run: cd /tmp/poore && python3 /path/to/src/prep.py  (needs voom.csv, meta.csv there)."""
import pandas as pd,numpy as np
ch=[];idx=[]
for c in pd.read_csv('voom.csv',index_col=0,chunksize=2000):
    idx+=list(c.index); ch.append(c.values.astype(np.float32)); cols=list(c.columns)
X=np.vstack(ch); print(X.shape); np.save('X.npy',X); pd.Series(idx).to_csv('idx.csv',index=False); pd.Series(cols).to_csv('cols.csv',index=False)
