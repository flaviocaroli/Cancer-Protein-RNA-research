# Headless version of notebook 02 for the smoke test (synthetic fallback data).
import sys, json
from pathlib import Path
import numpy as np, pandas as pd, joblib
sys.path.insert(0, str(Path.cwd()/"src"))
import harmonize as H, normalization as N

TRAIN=["CCLE","TMT","iTRAQ"]; TEST=["LabelFree"]; ALL=TRAIN+TEST
HARM="quantile_normal"; MODEL_NORM="zscore"
SPLIT="_".join(d.lower() for d in TRAIN)+"__"+"_".join(d.lower() for d in TEST)
SD=Path("data/splits")/SPLIT; SD.mkdir(parents=True, exist_ok=True)

def synth(name):
    rng=np.random.default_rng(abs(hash(name))%2**32)
    n={"CCLE":120,"TMT":60,"iTRAQ":40,"LabelFree":45}[name]
    genes=[f"P{i:03d}" for i in range(20)]+[f"G{i:04d}" for i in range(180)]; prots=[f"P{i:03d}" for i in range(20)]
    sc={"CCLE":1,"TMT":1.4,"iTRAQ":.7,"LabelFree":3}[name]; sh={"CCLE":0,"TMT":.5,"iTRAQ":-.4,"LabelFree":2}[name]
    W={p:rng.normal(size=6) for p in prots}; drv={p:rng.choice(len(genes),6,False) for p in prots}
    X=pd.DataFrame(rng.normal(sh,sc,(n,len(genes))),index=[f"{name}::s{i}" for i in range(n)],columns=genes)
    Y=pd.DataFrame(index=X.index,columns=prots,dtype=float)
    for p in prots: Y[p]=X.iloc[:,drv[p]].to_numpy()@W[p]+rng.normal(scale=.5,size=n)
    return X,Y
raw={d:dict(zip(["RNA","Protein"],synth(d))) for d in ALL}
cg=sorted(set.intersection(*[set(raw[d]["RNA"].columns) for d in ALL]))
cp=sorted(set.intersection(*[set(raw[d]["Protein"].columns) for d in ALL]))
tr_rna=pd.concat([raw[d]["RNA"][cg] for d in TRAIN]); tr_prot=pd.concat([raw[d]["Protein"][cp] for d in TRAIN])
rp=H.fit_reference_normalizer(tr_rna,method=HARM); pp=H.fit_reference_normalizer(tr_prot,method=HARM,transform="asinh")
harm={d:{"RNA":H.apply_reference_normalizer(raw[d]["RNA"][cg],rp),
         "Protein":H.apply_reference_normalizer(raw[d]["Protein"][cp],pp)} for d in ALL}
fn=N.fit_feature_normalizer(pd.concat([harm[d]["RNA"] for d in TRAIN])[cg],method=MODEL_NORM)
ap=lambda dfX: pd.DataFrame(N.apply_feature_normalizer(dfX[cg],fn),index=dfX.index,columns=cg)
Xtr=ap(pd.concat([harm[d]["RNA"] for d in TRAIN])); Ytr=pd.concat([harm[d]["Protein"] for d in TRAIN])[cp]
joblib.dump(Xtr,SD/"X_train.joblib"); joblib.dump(Ytr,SD/"Y_train.joblib")
for d in TEST:
    joblib.dump(ap(harm[d]["RNA"]),SD/f"X_test__{d}.joblib"); joblib.dump(harm[d]["Protein"][cp],SD/f"Y_test__{d}.joblib")
joblib.dump(cg,SD/"common_genes.joblib"); joblib.dump(cp,SD/"common_proteins.joblib")
joblib.dump({"rna":rp,"protein":pp,"feature_norm":fn},SD/"harmonization_params.joblib")
(SD/"split_config.json").write_text(json.dumps({"split_name":SPLIT,"train":TRAIN,"test":TEST,
  "harmonization":HARM,"model_norm":MODEL_NORM,"n_genes":len(cg),"n_proteins":len(cp)},indent=2))
print("built split", SPLIT, "X_train",Xtr.shape,"proteins",len(cp))
