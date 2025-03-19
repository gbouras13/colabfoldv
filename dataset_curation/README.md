# Dataset Construction

## RefSeq

* Downloaded all RefSeq virus proteins as of 13-09-24 from https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Protein&SourceDB_s=RefSeq 
* Total of `680608` proteins

## PhageScope


* Total of `873717` phage genomes were taken from [PhageScope](https://doi.org/10.1093/nar/gkad979) (under the heading 'Phage FASTA File Download' [here](https://phagescope.deepomics.org/download))
* Once combined, they were chunked into 999 batches with seqki


```bash
./seqkit split2 -p 999 all_phage_scope.fasta
```

* Then run Pharokka v1.7.3 was run with `pyrodigal-gv` on each of the 999 batches using the [Biocontainers singularity container](http://www.biocontainers.com/pharokka?tab=tags&tag=latest)

```
 singularity exec  $containerImage \
         pharokka.py -d pharokka_db -t $THREADS -g prodigal-gv -m -i all_phage_scope.fasta.split/all_phage_scope.part_${formatted_input}.fasta  -o $OUTDIR/all_phage_scope.part_${formatted_input} -f --skip_extra_annotations --skip_mash  
```


* Once combined, a total of `43632192` proteins were called in these genomes

## enVhogs

* The [enVhogs](https://www.biorxiv.org/content/10.1101/2024.06.25.600602v2) are viral protein families from metagenomes. The "new" enVhogs can be accessed from http://envhog.u-ga.fr/envhog/. However, I did the following on 31 August 2024 with the "old" version (that have now been removed and updated) - if you are desperate to reproduce this for some reason, make an issue. 
* This in hindsight was done pretty inefficiently so bear with me
* Overall, I took `25543277` proteins from the "old" enVhogs as follows:

* With `EnVhog.tar.gz` untar and unpack the MSAs

```bash
tar -xzf EnVhog.tar.gz

# 3. convert the hhsuite format to get a3m files
# install hhsuite v3.3
conda create -n hhsuite hhsuite
conda activate hhsuite

# get a3m files using ffindex_unpack
mkdir -p a3ms
ffindex_unpack envhog_hmm/EnVhog_a3m.ffdata  envhog_hmm/EnVhog_a3m.ffindex  a3ms/ .
```

* Then run `python put_subdirectory.py` to distribute the MSAs into subdirectories, `python filter_fastas.py` to extract the protein sequence from member of each MSA (not consensus sequences) and finally `concat_fastas.py` to combine them all into one file

```bash
python put_subdirectory.py

```

## Jaeger

* Phages and prophage were predicted from Mgnify contigs using [Jaeger](https://github.com/Yasas1994/Jaeger) and shared with me. These were then annotated with Pharokka v1.7.3 and gene called with `pyrodigal-gv` in the same way as `Phage scope`
* In total there were `174561888` proteins from this source

## Deduplication

* All proteins were combined, yielding `244417965` proteins ( `combined_virus_proteins.faa` )
* MMseqs2 v `87e7103d289029dc3345f85ea9a4c4c6d6416e46` was then used to deduplicate, cluster (at 30% seq-id and 90% coverage) and create colabfold compatible databases as follows
* Specifically, the database contains `129944764` non-redundant proteins and `60287451` clusters

```
THREADS=256
mmseqs createdb combined_virus_proteins.faa  combined_virus_proteins_DB
mmseqs linclust combined_virus_proteins_DB  combined_virus_proteins_DB100 tmp --min-seq-id 1.0  -c 1.0 --threads $THREADS
mmseqs createclusearchdb  combined_virus_proteins_DB  combined_virus_proteins_DB100  combined_virus_proteins_DB100_search  --threads $THREADS
mmseqs linclust combined_virus_proteins_DB100_search  combined_virus_proteins_DB100_clu tmp --min-seq-id 0.3  -c 0.9 --threads $THREADS
mmseqs createclusearchdb  combined_virus_proteins_DB100_search  combined_virus_proteins_DB100_clu  combined_virus_proteins_DB100_clu_searchDB  --threads $THREADS
mmseqs align combined_virus_proteins_DB100_search  combined_virus_proteins_DB100_search combined_virus_proteins_DB100_clu  combined_virus_proteins_DB100_aln -e inf -a --threads $THREADS
```