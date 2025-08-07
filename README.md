# colabfoldv

This repository is intended to replicate only the `colabfold_search` functionality of the local ColabFold MSA generation pipeline, with an added protein database of `129944764` non-redundant proteins tailored for viruses, especially phages. Please see `dataset_curation` for more details on how it was constructed.

For the full Colabfold functionality and up-to-date changes (and frankly anything outside of the use case presented below), please go to the [ColabFold repository](https://github.com/sokrypton/ColabFold).

This was used to create some MSAs used in [phold's](https://github.com/gbouras13/phold) database.

## Installation

* This will replicate the MSAs used in the enVhogs and PHROG singleton protein structures in [Phold](https://github.com/gbouras13/phold). You can use other, later versions of MMSeqs2 if you'd like.

* This assumes you have [conda](https://github.com/conda-forge/miniforge) installed.

```bash
conda create -n colabfoldv_env python=3.12 mmseqs==15.6f452 pip
conda activate colabfoldv_env
git clone https://github.com/gbouras13/colabfoldv
cd colabfoldv
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

### Download the viral/phage DB

* The viral/phage database database is stored on Zenodo. It is 39GB, and therefore it may take some time to download.
* I would highly recommend using aria2c (as used in MMSeqs2 and Foldseek) to download, but you can use e.g. wget or curl

```bash
aria2c https://zenodo.org/records/15045387/files/viral_db.tar.xz
tar -xvf viral_db.tar.xz
```

### To use (change to the desired number of threads)

* This will create MSAs without pairing (i.e. for monomers) for all proteins in `example/NC_043029_aa.fasta` using all three databases

```bash
THREADS=72
colabfold_search example/NC_043029_aa.fasta colabfoldDBs viral_db NC_043029_aa_msas --threads $THREADS
```

