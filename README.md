# Heart
Boolean Analysis of Heart Datasets

## CST-KO Heart Analysis

Title: Catestatin improves cardiac metabolic flexibility by reversing insulin
resistance and potentiates mitochondrial function

To reproduce the results:
1. Clone the github repositories github:sahoo00/Hegemon and sahoo00/BoNE
2. Download the data
[GSE119087](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119087) and
[BoLHD Dataset](https://zenodo.org/records/19722879).

3. Organize the data folder
```bash
$ ls data/GSE119087_human-gpl570-*
data/GSE119087_human-gpl570-bv.txt    data/GSE119087_human-gpl570-ih.txt
data/GSE119087_human-gpl570-idx.txt   data/GSE119087_human-gpl570-thr.txt
data/GSE119087_human-gpl570-expr.txt
```
4. Run the Notebook [CST-KO-Heart.ipynb](CST-KO-Heart.ipynb)

### BoLHD Dataset

We created Boolean Lab Heart Disease Dataset (BoLHD Dataset) and submitted in
the following repository.
Download this dataset from public repository[1].
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19722879.svg)](https://doi.org/10.5281/zenodo.19722879)


### BoLHD Analysis

```bash
tar xvzf BoLHD.tar.gz
mv dataset data
```

Please refer to the jupyter notebook [CST-KO-Heart.ipynb](CST-KO-Heart.ipynb) to reproduce the results.

## References
[1] Sahoo, D. (2026). Boolean Lab Heart Disease Dataset (Version v1) [Data set].
Zenodo. https://doi.org/10.5281/zenodo.19722879
