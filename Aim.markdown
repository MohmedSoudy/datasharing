# Analysis to Extract the transposable elements from WGS reads

## TEMP2
```
sudo apt update 
sudo apt-get update
#Install bwa
sudo apt-get -y install bwa
#Install samtools
#Install HTSLIB
git clone https://github.com/samtools/htslib.git
autoreconf -i  # Build the configure script and install files it uses
./configure    # Optional but recommended, for choosing extra functionality
make
make install
sudo apt-get update
sudo apt-get install gcc
sudo apt-get install make
sudo apt-get install libbz2-dev
sudo apt-get install zlib1g-dev
sudo apt-get install libncurses5-dev 
sudo apt-get install libncursesw5-dev
sudo apt-get install liblzma-dev
#Install samtools
cd ..
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar -vxjf samtools-1.9.tar.bz2
cd samtools-1.9
make
#Install BCFTools
cd ..
wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2
tar -vxjf bcftools-1.9.tar.bz2
cd bcftools-1.9
make
#Install bedtobigbed
#Install Anaconda
curl –O https://repo.anaconda.com/archive/Anaconda3-2022.05-Linux-x86_64.sh
sha256sum Anaconda3–2020.02–Linux–x86_64.sh
bash Anaconda3-2022.05-Linux-x86_64.sh
conda install -c bioconda ucsc-bedtobigbed
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
tar -zxvf bedtools-2.29.1.tar.gz
cd bedtools2
make
wget 
#Install trinity
sudo apt-get update
sudo apt-get install trinity
source ~/.bashrc
#Install TEMP2
git clone https://github.com/weng-lab/TEMP2
ln -s $PWD/TEMP2/TEMP2 /usr/bin
#Run TEMP2
cd TEMP2/test
TEMP2 insertion -l test.1.fastq.gz -r test.2.fastq.gz -I bwa_index/chr2L -g chr2L.fa -R transposon.fa -t rmsk.bed -o test_output -c 2
#Output is listed in the TEMP2 Output folder as follows 
** test.insertion.bed: a bed format file with additional 8 columns. **
- Column 1,2&3: Reference genome position (chromosome, start, and end) of the transposon insertion.
- Column 4: Information of inserted transposon, including transposon name, start coordinate and end coordinate of the inserted transposon, and which strand it inserts. Separated by :
- Column 5: Frequency of the inserted transposon. It generally means what fraction of sequenced genome present this insertion.
- Column 6: Which strand this transposon inserts.
- Column 7: Type of the insertion. We typically separate insertions into three categories based on how many reads support them.
- Column 8: Number of reads supporting this insertion.
- Column 9: Number of reads that do not support this insertion, AKA reference reads.
- Column 10: Number of supporting reads at 5'end of the insertion.
- Column 11: Number of supporting reads at 3'end of the insertion.
- Column 12: Target site duplication (TSD) of the insertion. unknown is shown if not applicable.
- Column 13: Reliability of this insertion (0–100). 100 for 2p and 1p1 insertions. For singleton insertions, TEMP2 already filtered out most of the false positives but not all of them. The reliability is a percentage stand for how many singleton insertions of a specific transposon is.
- Column 14: Number of supporting reads at 5'end of the insertion junction.
- Column 15: Number of supporting reads at 3'end of the insertion junction.
** test.soma.summary.txt: a tab delimited file includes 10 columns ** 
- Column 1: Name of transposon.
- Column 2: Estimated number of de novo insertion of this transposon per genome.
- Column 3: 95th percentile (lambda distribution) number of de novo insertion of this transposon per genome.
- Column 4: Total estimated number of de novo insertion of this transposon.
- Column 5: 95th percentile (lambda distribution) number of de novo insertion of this transposon.
- Column 6: Number of singleton reads mapped to the end regions (positive regions) of this transposon.
- Column 7: Number of singleton reads mapped to the center region (negative region) of this transposon.
- Column 8: Number of reads mapped to the end regions (positive regions) of this transposon.
- Column 9: Number of reads mapped to the center region (negative region) of this transposon.
- Column 10: Status of the estimation.
** test.supportReadsUnfiltered.bb: bigBed6 file for all supporting reads.**
** test.supportingRead.dis.pdf: figures showing where the supporting reads mapped to each transposon.**
 ```
## InMut-Finder 

```
git clone https://github.com/jsg200830/InMut-Finder.git
./run_command.sh
#Output is listed in the InMut-Finder_Output folder 
- The first one contains the summary of insertion sites, with five columns, which indicate the information of location (genomic coordinates of insertion), gap (extend of insertion), counts (read counts to support this insertion), gene (neighbor genes within a distance of 2000 bp) and P-value (probability value to reject the hypothesis in Poisson distribution). See file of “test.fa_insertion_list_nu.txt_score.csv”.

- The second contains the 5’ and 3’ sequences from the Nanopore/PacBio sequencing for each insertion. The sequences could be used for designing primers to validate the insertion. See file of “test.fa_insertion_flanking_seqs_nu.txt”.

- The third contains all the reads which cover both the inserted target fragment and the flanking sequences which can be mapped to the reference genome. This file is useful to do the alignment against the reference genome, for example by minimap2. The bam file can be viewed in IGV tool for manually examining the insertion sites. See file of “test.fa_flanking.fasta”.
```
## xTea

```
git clone https://github.com/parklab/xTea.git    #latest version
#pre-processed repeat library used by xTea (this library is used for both short and long reads)
wget https://github.com/parklab/xTea/raw/master/rep_lib_annotation.tar.gz
#Download the gff enome annotation from https://www.gencodegenes.org/human/release_33.html
#Install pysam
conda config --add channels r
conda config --add channels bioconda
conda install pysam -y
conda install sortedcontainers -y    #Install sortedcontainers
python -m pip install deep-forest
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
#install xTea
conda install -y xtea
#replace sample id in sample.txt with our sample id (NA24385)
xtea -i sample_id.txt -b illumina_bam_list.txt -x null -p ./path_work_folder/ -o submit_jobs.sh -l /home/rep_lib_annotation/ -r /home/reference/genome.fa -g /home/gene_annotation_file.gff3 --xtea /mnt/e/Testing-Project/xTea/xtea/ -f 5907 -y 7  --slurm -t 0-12:00 -q short -n 8 -m 25
./run_gnrt_pipeline.sh 
./submit_jobs.sh 
```
