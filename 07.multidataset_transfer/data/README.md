# data/
Notebook 00 writes per-platform matrices to `loaded/`:
  <D>_rna.joblib, <D>_protein.joblib for D in {CCLE,TMT,iTRAQ,LabelFree}
Notebook 02 writes each split to `splits/<split_name>/` (matrices + params + split_config.json).
Place Collection A under data/Collection_A/ and CCLE joblibs as data/CCLE_rna.joblib / CCLE_protein.joblib.
No real data? The notebooks fall back to synthetic; `run_smoke_test.sh` builds a synthetic split.
