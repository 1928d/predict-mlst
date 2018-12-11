# predict-mlst
In this exercise, we download the allelic sequences for the 7 MLST genes for E. coli,
from the pubMLST database. Then, we download the E. coli reference genome, from NCBI genome database.
We blast all the gene alleles to the reference genome, and we download the mlst profiles
for E. coli from pubMLST. The mlst profiles are the list of strain types and the allelic combinations
of the 7 genes that correspond to them.
Finally we use a script that parses the blast result and the file with the mlst profiles
and we get the strain type of the reference genome.

## if you don't have blast installed, do as following:
```
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz
tar xzf ncbi-blast-2.6.0+-x64-linux.tar.gz
nano .bashrc
```
scroll to the end of the file and add:
```
export PATH="$PATH:/home/alex/ncbi-blast-2.6.0+/bin"
```
then, after quiting the nano editor:
```
source .bashrc
```

## create analysis folder
```
mkdir mlst_analysis
cd mlst_analysis
```

## get mlst gene alleles
```
wget https://pubmlst.org/data/alleles/ecoli/{adk,fumC,gyrB,icd,mdh,purA,recA}.tfa
```

## tidy up the folder
```
cat *tfa > mlst_alleles.fna
rm *tfa
```

## get ecoli assembled genome
```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz
```

## convert it into a blast parseable database
```
makeblastdb -in GCF_000005845.2_ASM584v2_genomic.fna -dbtype nucl -out ecoli
```

## blast mlst gene alleles to ecoli genome
```
blastn -query mlst_alleles.fna -db ecoli -out mlst_result.txt -outfmt 6 -qcov_hsp_perc 100 -perc_identity 100
```

## get mlst profiles
```
wget https://pubmlst.org/data/profiles/ecoli.txt
```

## use the script to get the strain type
```
python3 predict_mlst_profile.py -b mlst_result.txt -p ecoli.txt
```
