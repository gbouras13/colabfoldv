# colabfoldv

This repository is intended to replicate only the `colabfold_search` functionality of the local ColabFold MSA generation pipeline, with an added protein database of `129944764` proteins tailored for viruses, especially phages. Please see `dataset_curation` for more details on how it was constructed.

For the full Colabfold functionality and up-to-date changes (and frankly anything outside of the use case presented below), please go to the [ColabFold repository](https://github.com/sokrypton/ColabFold).

## Installation

* This will replicate the MSAs used in the enVhogs and PHROG singleton protein structures in [Phold](https://github.com/gbouras13/phold). You can use other, later versions of MMSeqs2 if you'd like.

```bash
conda create -n colabfoldv_test python=3.12 mmseqs==15.6f452
conda activate colabfoldv_test
pip install -e .
pip install colabfold[alphafold]
```

### Download Colabfold DBs

* This will download the regular ColabFold uniref30_2302 and colabfold_envdb_202108 databases

```bash
mkdir -p colabfoldDBs
cd colabfoldDBs
bash ../setup_databases.sh
```

### Download the viral DB

* The viral database database is stored on Zenodo (just). It therefore may take some time to download.

```bash

```

### To use (change to the desired number of threads)


```bash
THREADS=72
colabfold_search example/NC_043029_aa.fasta colabfoldDBs viral_db NC_043029_aa_msas --threads $THREADS
```

