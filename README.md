# This is my workflow in Linux sever

## Working Environment 

```
[wirawit.s@omic Y_Lipolytica]$ tree -L 1
.
├── 01.raw_data
├── 02.QC
├── 03.trimming
├── 04.QC_after
├── 05.mapping_reads
├── 06.quantify
├── _archived
├── multi_qc_all.sh
└── PBS-script.sh

7 directories, 2 files
```

## Directory tree
```
[wirawit.s@omic Y_Lipolytica]$ tree -L 2
.
├── 01.raw_data
│   ├── ESIG
│   ├── NOV24TCS03076_01_X401SC24032711-Z01-F001_C
│   ├── NOV24TCS03076_01_X401SC24032711-Z01-F001_C.zip
│   └── NOV24TCS03076_01_X401SC24032711-Z01-F001_C.zipZone.Identifier
├── 02.QC
│   ├── CXAUAI_qc
│   └── ESIG_QC
├── 03.trimming
│   ├── compressing.sh
│   ├── cutadapttest
│   ├── CXAUAI
│   ├── ESIG
│   ├── NOV24TCS03076_01_X401SC24032711-Z01-F001_C
│   ├── parallel.err
│   ├── parallel.log
│   ├── Ref_Genome
│   └── test_hisat.sh
├── 04.QC_after
│   ├── CXAUAI_QC_after
│   ├── ESIG_QC_after
│   ├── multiqc_data
│   └── multiqc_report.html
├── 05.mapping_reads
│   ├── CXAUAI
│   ├── ESIG
│   ├── genome_index
│   ├── genome_index_notrim
│   ├── multiqc_result
│   ├── Ref_genome
│   └── Test
├── _archived
│   ├── ESIG
│   ├── NOV24TCS03076_01_X401SC24032711-Z01-F001_C
│   ├── NOV24TCS03076_01_X401SC24032711-Z01-F001_C.zip
│   └── NOV24TCS03076_01_X401SC24032711-Z01-F001_C.zipZone.Identifier
├── PBS-script.sh
└── workflow.txt

27 directories, 11 files
```

### Note
Don't be bothered by files with ?Zone.Identifier, it happened when I moved files from Windows to WSL.

## 1. Upload files

### 1.1 Upload files with FileZilla from my WSL

For ESIG (mutant):
```
[wirawit.s@omic Y_Lipolytica]$ tree -L 1 _archived/ESIG/
_archived/ESIG/
├── ESIG_SE-1
├── ESIG-SE-2
├── ESIG-SE-3
├── ESIG-SG-1
├── ESIG-SG-2
├── ESIG-SG-3
└── Reference genome

7 directories, 0 files
```

For CXAUAI (wild-type):
```
[wirawit.s@omic Y_Lipolytica]$ tree -L 1 _archived/ESIG/
_archived/ESIG/
├── ESIG_SE-1
├── ESIG-SE-2
├── ESIG-SE-3
├── ESIG-SG-1
├── ESIG-SG-2
├── ESIG-SG-3
└── Reference genome

7 directories, 0 files
```
```
[wirawit.s@omic Y_Lipolytica]$ tree -L 2 _archived/NOV24TCS03076_01_X401SC24032711-Z01-F001_C
_archived/NOV24TCS03076_01_X401SC24032711-Z01-F001_C
├── 01.RawData
│   ├── CXAUAI_SE1
│   ├── CXAUAI_SE2
│   ├── CXAUAI_SE3
│   ├── CXAUAI_SG1
│   ├── CXAUAI_SG2
│   └── CXAUAI_SG3
├── 02.Report_X401SC24032711-Z01-F001.zip
├── checkSize.xls
├── MD5-manual.pdf
├── MD5.txt
├── MD5-win.exe
├── Readme.html
└── X401SC24032711-Z01-F001.tar

7 directories, 7 files
```

### 1.2 Run md5checksum

For ESIG (mutant):
Example: go to each directory
```
[wirawit.s@omic ESIG]$ cd ESIG_SE-2/
[wirawit.s@omic ESIG_SE-2]$ ls *
ESIG-SE-2_1.fq.gz                  ESIG-SE-2_2.fq.gz?Zone.Identifier
ESIG-SE-2_1.fq.gz?Zone.Identifier  MD5.txt
ESIG-SE-2_2.fq.gz                  MD5.txt?Zone.Identifier
```
Run `md5sum -c MD5.txt` in every directory manually
```
[wirawit.s@omic ESIG_SE-2]$ md5sum -c MD5.txt
ESIG-SE-2_1.fq.gz: OK
ESIG-SE-2_2.fq.gz: OK
```

For CXAUAI (wild-type):
Change directory to `NOV24TCS03076_01_X401SC24032711-Z01-F001_C/` (decompress first if needed)
```
[wirawit.s@omic 01.raw_data]$ cd NOV24TCS03076_01_X401SC24032711-Z01-F001_C/
[wirawit.s@omic NOV24TCS03076_01_X401SC24032711-Z01-F001_C]$ ls
02.Report_X401SC24032711-Z01-F001.zip  MD5-manual.pdf  MD5-win.exe
checkSize.xls                          MD5new.txt      Readme.html
CXAUAI                                 MD5.txt         X401SC24032711-Z01-F001.tar
```
This is more convenient than ESIG, but there's a catch: since I changed the `01.RawData` directory containing CXAUAI pairs to `CXAUAI`, I have to change the directory in `MD5.txt` to a new `MD5new.txt`.
```
[wirawit.s@omic NOV24TCS03076_01_X401SC24032711-Z01-F001_C]$ md5sum -c MD5new.txt
CXAUAI/CXAUAI_SE2/CXAUAI_SE2_1.fq.gz: OK
CXAUAI/CXAUAI_SE2/CXAUAI_SE2_2.fq.gz: OK
CXAUAI/CXAUAI_SE3/CXAUAI_SE3_1.fq.gz: OK
CXAUAI/CXAUAI_SE3/CXAUAI_SE3_2.fq.gz: OK
CXAUAI/CXAUAI_SE1/CXAUAI_SE1_1.fq.gz: OK
CXAUAI/CXAUAI_SE1/CXAUAI_SE1_2.fq.gz: OK
CXAUAI/CXAUAI_SG1/CXAUAI_SG1_1.fq.gz: OK
CXAUAI/CXAUAI_SG1/CXAUAI_SG1_2.fq.gz: OK
CXAUAI/CXAUAI_SG2/CXAUAI_SG2_1.fq.gz: OK
CXAUAI/CXAUAI_SG2/CXAUAI_SG2_2.fq.gz: OK
CXAUAI/CXAUAI_SG3/CXAUAI_SG3_1.fq.gz: OK
CXAUAI/CXAUAI_SG3/CXAUAI_SG3_2.fq.gz: OK
02.Report_X401SC24032711-Z01-F001.zip: OK
```

## 2. Run fastqc and multiqc

```
qsub -I -l nodes=1:ppn=1 -l walltime=10:00:00 -l mem=10gb #(not sure on ppn= and mem=)
module load bio/fastqc/0.11.9
module load bio/MultiQC/1.9-foss-2020a-Python-3.8.2
```

### Go to `
~/Y_Lipolytica/01.raw_data/ESIG`

```
cd ~/Y_Lipolytica/01.raw_data/ESIG
[wirawit.s@omic ESIG]$ ls
ESIG_SE-1  ESIG_SE-2  ESIG_SE-3  ESIG_SG-1  ESIG_SG-2  ESIG_SG-3
```

Run fastqc and specify output to `~/Y_Lipolytica/02.QC/ESIG_QC/ALL_sametime`
```
fastqc -threads 12 ESIG _S*/*.fq.gz -o ~/Y_Lipolytica/02.QC/ESIG_QC/ALL_sametime
```
Change directory to `~/Y_Lipolytica/02.QC/ESIG_QC/ALL_sametime `and run
```
multiqc *_fastqc.zip
```

here are the result
```
[wirawit.s@omic ALL_sametime]$ ls
ESIG-SE-1_1_fastqc.html  ESIG-SE-3_1_fastqc.zip   ESIG-SG-2_2_fastqc.html
ESIG-SE-1_1_fastqc.zip   ESIG-SE-3_2_fastqc.html  ESIG-SG-2_2_fastqc.zip
ESIG-SE-1_2_fastqc.html  ESIG-SE-3_2_fastqc.zip   ESIG-SG-3_1_fastqc.html
ESIG-SE-1_2_fastqc.zip   ESIG-SG-1_1_fastqc.html  ESIG-SG-3_1_fastqc.zip
ESIG-SE-2_1_fastqc.html  ESIG-SG-1_1_fastqc.zip   ESIG-SG-3_2_fastqc.html
ESIG-SE-2_1_fastqc.zip   ESIG-SG-1_2_fastqc.html  ESIG-SG-3_2_fastqc.zip
ESIG-SE-2_2_fastqc.html  ESIG-SG-1_2_fastqc.zip   multiqc_data
ESIG-SE-2_2_fastqc.zip   ESIG-SG-2_1_fastqc.html  multiqc_report.html
ESIG-SE-3_1_fastqc.html  ESIG-SG-2_1_fastqc.zip
```

### Go to
`~/Y_Lipolytica/01.raw_data/NOV24TCS03076_01_X401SC24032711-Z01-F001_C/CXAUAI`

Run fastqc and specified output to `~/Y_Lipolytica/02.QC/CXAUAI_qc/ALL`

```
fastqc -threads 12 CXAUAI_S*/*.fq.gz -o ~/Y_Lipolytica/02.QC/CXAUAI_qc/ALL
```

& change directory to `~/Y_Lipolytica/02.QC/CXAUAI_qc/ALL` and run

```
multiqc *_fastqc.zip
```
here are the result likewise ESIG
```
[wirawit.s@omic ALL]$ ls
CXAUAI_SE1_1_fastqc.html  CXAUAI_SE3_1_fastqc.zip   CXAUAI_SG2_2_fastqc.html
CXAUAI_SE1_1_fastqc.zip   CXAUAI_SE3_2_fastqc.html  CXAUAI_SG2_2_fastqc.zip
CXAUAI_SE1_2_fastqc.html  CXAUAI_SE3_2_fastqc.zip   CXAUAI_SG3_1_fastqc.html
CXAUAI_SE1_2_fastqc.zip   CXAUAI_SG1_1_fastqc.html  CXAUAI_SG3_1_fastqc.zip
CXAUAI_SE2_1_fastqc.html  CXAUAI_SG1_1_fastqc.zip   CXAUAI_SG3_2_fastqc.html
CXAUAI_SE2_1_fastqc.zip   CXAUAI_SG1_2_fastqc.html  CXAUAI_SG3_2_fastqc.zip
CXAUAI_SE2_2_fastqc.html  CXAUAI_SG1_2_fastqc.zip   multiqc_data
CXAUAI_SE2_2_fastqc.zip   CXAUAI_SG2_1_fastqc.html  multiqc_report.html
CXAUAI_SE3_1_fastqc.html  CXAUAI_SG2_1_fastqc.zip
```


## 3.Trim_Adapter with Cutadapt

### Change Directory to `~/Y_Lipolytica/03.trimming/`

I copied `ESIG` and `NOV24TCS03076_01_X401SC24032711-Z01-F001_C` directories included files inside from `Y_Lipolytica/01.raw_data`
with
```
cp -r ../01.raw_data/ESIG/ .
```
and
```
cp -r ../01.raw_data/NOV24TCS03076_01_X401SC24032711-Z01-F001_C/ .
```
then move ` CXAUAI` in `NOV24TCS03076_01_X401SC24032711-Z01-F001_C` to current directory

directory list :
```
[wirawit.s@omic 03.trimming]$ ls
CXAUAI  NOV24TCS03076_01_X401SC24032711-Z01-F001_C ESIG 
```
Start with Wild-type(s)
these are file within `CXAUAI` at first you should't had .sh,.err,.log and files begin with "trimmed"

```
[wirawit.s@omic 03.trimming]$ tree CXAUAI/
CXAUAI/
├── CXAUAI_SE1
│   ├── CXAUAI_SE1_1.fq.gz
│   ├── CXAUAI_SE1_2.fq.gz
│   ├── MD5.txt
│   ├── trimmed_CXAUAI_SE1_1.fq.gz
│   └── trimmed_CXAUAI_SE1_2.fq.gz
├── CXAUAI_SE2
│   ├── CXAUAI_SE2_1.fq.gz
│   ├── CXAUAI_SE2_2.fq.gz
│   ├── MD5.txt
│   ├── trimmed_CXAUAI_SE2_1.fq.gz
│   └── trimmed_CXAUAI_SE2_2.fq.gz
├── CXAUAI_SE3
│   ├── CXAUAI_SE3_1.fq.gz
│   ├── CXAUAI_SE3_2.fq.gz
│   ├── MD5.txt
│   ├── trimmed_CXAUAI_SE3_1.fq.gz
│   └── trimmed_CXAUAI_SE3_2.fq.gz
├── CXAUAI_SG1
│   ├── CXAUAI_SG1_1.fq.gz
│   ├── CXAUAI_SG1_2.fq.gz
│   ├── MD5.txt
│   ├── trimmed_CXAUAI_SG1_1.fq.gz
│   └── trimmed_CXAUAI_SG1_2.fq.gz
├── CXAUAI_SG2
│   ├── CXAUAI_SG2_1.fq.gz
│   ├── CXAUAI_SG2_2.fq.gz
│   ├── MD5.txt
│   ├── trimmed_CXAUAI_SG2_1.fq.gz
│   └── trimmed_CXAUAI_SG2_2.fq.gz
├── CXAUAI_SG3
│   ├── CXAUAI_SG3_1.fq.gz
│   ├── CXAUAI_SG3_2.fq.gz
│   ├── MD5.txt
│   ├── trimmed_CXAUAI_SG3_1.fq.gz
│   └── trimmed_CXAUAI_SG3_2.fq.gz
├── CXAUAI_trimming.sh
├── parallel.err
└── parallel.log

6 directories, 33 files

```
Then i submit `CXAUAI_trimming.sh` with following command
```
qsub CXAUAI_trimming.sh 
```
here is detail of the script
important part is .sh script has to be in  the same directory with CXAUAI_S... directory (`~/Y_Lipolytica/03.trimming/CXAUAI`)

```
### Declare job non-rerunable
#PBS -r n
### Output files
#PBS -e parallel.err
#PBS -o parallel.log
### Mail to user
#PBS -m ae
#PBS -l walltime=96:00:00
#PBS -l nodes=1:ppn=24
#PBS -l mem=100gb
#################################################
#This job's working directory
echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR
###############################################

==THIS PART required for submit jobs with qsub command==

```
==THIS PART==   required for submit jobs with  `qsub` command located in most of scripts in I write 

```sh

#!/bin/bash
module load bio/cutadapt/2.10-GCCcore-9.3.0-Python-3.8.2
# Iterate over directories
for dir in */; do
    # Get the directory name
    NAMES="${dir%/}"

    # Define the command to be executed
    COMMAND="cutadapt -j 0 -m 60 -u 15 -U 15 -a AGATCGGAAGAG -A AGATCGGAAGAG -o trimmed_${NAMES}_1.fq.gz -p trimmed_${NAMES}_2.fq.gz *_1.fq.gz *_2.fq.gz"

    echo "Running '$COMMAND' in directory: $dir"
    (cd "$dir" && eval "$COMMAND")
done


###############################################

```

and likewise for Mutant .sh script has to be in  the same directory with ESIG_S... dir (`~/Y_Lipolytica/03.trimming/ESIG`)

```
[wirawit.s@omic 03.trimming]$ tree ESIG/
ESIG/
├── ESIG_SE-1
│   ├── ESIG-SE-1_1.fq.gz
│   ├── ESIG-SE-1_2.fq.gz
│   ├── MD5.txt
│   ├── trimmed_ESIG_SE-1_1.fq.gz
│   └── trimmed_ESIG_SE-1_2.fq.gz
├── ESIG_SE-2
│   ├── ESIG-SE-2_1.fq.gz
│   ├── ESIG-SE-2_2.fq.gz
│   ├── MD5.txt
│   ├── trimmed_ESIG_SE-2_1.fq.gz
│   └── trimmed_ESIG_SE-2_2.fq.gz
├── ESIG_SE-3
│   ├── ESIG-SE-3_1.fq.gz
│   ├── ESIG-SE-3_2.fq.gz
│   ├── MD5.txt
│   ├── trimmed_ESIG_SE-3_1.fq.gz
│   └── trimmed_ESIG_SE-3_2.fq.gz
├── ESIG_SG-1
│   ├── ESIG-SG-1_1.fq.gz
│   ├── ESIG-SG-1_2.fq.gz
│   ├── MD5.txt
│   ├── trimmed_ESIG_SG-1_1.fq.gz
│   └── trimmed_ESIG_SG-1_2.fq.gz
├── ESIG_SG-2
│   ├── ESIG-SG-2_1.fq.gz
│   ├── ESIG-SG-2_2.fq.gz
│   ├── MD5.txt
│   ├── trimmed_ESIG_SG-2_1.fq.gz
│   └── trimmed_ESIG_SG-2_2.fq.gz
├── ESIG_SG-3
│   ├── ESIG-SG-3_1.fq.gz
│   ├── ESIG-SG-3_2.fq.gz
│   ├── MD5.txt
│   ├── trimmed_ESIG_SG-3_1.fq.gz
│   └── trimmed_ESIG_SG-3_2.fq.gz
├── ESIG_trimming.sh
├── parallel.err
└── parallel.log

6 directories, 33 files
```

submit job

```
qsub ESIG_trimming.sh
```

here i script for Mutant

```sh

#!/bin/bash
module load bio/cutadapt/2.10-GCCcore-9.3.0-Python-3.8.2
# Iterate over directories
for dir in */; do
    # Get the directory name
    NAMES="${dir%/}"

    # Define the command to be executed
    COMMAND="cutadapt -j 0 -m 60 -u 15 -U 15 -a AGATCGGAAGAG -A AGATCGGAAGAG -o trimmed_${NAMES}_1.fq.gz -p trimmed_${NAMES}_2.fq.gz *_1.fq.gz *_2.fq.gz"

    echo "Running '$COMMAND' in directory: $dir"
    (cd "$dir" && eval "$COMMAND")
done


###############################################
```
As tree command result above ,files begin with "trim" will be generated as result for cutadapt trimming process


## (additional) QC after trimming

### run Fastqc  after trimming process
change directory to `~/Y_Lipolytica/03.trimming/CXAUAI`
```
cd ~/Y_Lipolytica/03.trimming/CXAUAI
fastqc -threads 12 CXAUAI_S*/*.fq.gz -o ../../04.QC_after/CXAUAI_QC_after
```
change directory to `~/Y_Lipolytica/03.trimming/ESIG`
```
cd ~/Y_Lipolytica/03.trimming/ESIG
fastqc -threads 12 ESIG_S*/**.fq.gz -o ../../04.QC_after/ESIG_QC_after
```

### run multiqc
change directory to `~/04.QC_after/ESIG_QC_after`
run multiqc
```
multiqc . 
```

## 4. Align reads to genome with STAR
### Change directory to `Y_Lipolytica/05.mapping_reads`

```
[wirawit.s@omic 05.mapping_reads]$ ls
CXAUAI          newgenomeindex  notes.txt     ReadPerGen_sum
ESIG            newscript.sh    parallel.err  Ref_genome
multiqc_result  new_withcounts  parallel.log  Test

```

### 4.1 Creating genome index folder
for nontrim fq and do the same for trimed

```sh
#!/bin/bash
module load star/2.7.6a-yk2el5n

  STAR --runThreadN 10 \
  --runMode genomeGenerate \
  --genomeSAindexNbases 11 \
  --genomeDir ~/Y_Lipolytica/05.mapping_reads/newgenomeindex \
  --genomeFastaFiles ~/Y_Lipolytica/05.mapping_reads/Ref_genome/GCF_00
  --sjdbGTFfile ~/Y_Lipolytica/05.mapping_reads/Ref_genome/GCF_0000025
  --sjdbOverhang 99

```
### 4.2mapping reads to genome

run four script
```
qsub 'script_name.sh'
```
i.e.  CXAUAI_n_trim
the part shown above are same for script that use `qsub` command
What these 4 scripts essentialy do are mapping 12 fq each to genome index
i.e. notrim_Wild-type

```sh
#!/bin/bash

# Set the directory where your data is stored
data_dir=~/Y_Lipolytica/03.trimming/CXAUAI/CXAUAI_S*

module load star/2.7.6a-yk2el5n
# Loop through each directory in the data directory
for dir in $data_dir; do
    echo "STAR is working on files in directory: $dir"
    echo "directory name: ${dir##*/}"
    
    #create directory for output
    mkdir -p ~/Y_Lipolytica/05.mapping_reads/CXAUAI/CXAUAI_n_trim/${dir##*/}
    
    # Set the file patterns
    pattern_1="$dir/CXAUAI_S*_1.fq.gz"
    pattern_2="$dir/CXAUAI_S*_2.fq.gz"
    
    # List files matching the patterns
    echo "Files matching pattern 1:"
    echo $pattern_1
    
    echo "Files matching pattern 2:"
    echo $pattern_2


    STAR --runThreadN 24 \
    --genomeDir ~/Y_Lipolytica/05.mapping_reads/newgenomeindex  \
    --readFilesIn $pattern_1 $pattern_2 \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --outFileNamePrefix ~/Y_Lipolytica/05.mapping_reads/CXAUAI/CXAUAI_n_trim/${dir##*/}/${dir##*/} \
    --limitBAMsortRAM 50000000000

    
done 

```


### Result
i.e. CXAUAI_n_trim
```
[wirawit.s@omic CXAUAI_n_trim]$ tree -L 2
.
├── CXAUAI_n_trim.sh
├── CXAUAI_SE1
│   ├── CXAUAI_SE1Aligned.sortedByCoord.out.bam
│   ├── CXAUAI_SE1Log.final.out
│   ├── CXAUAI_SE1Log.out
│   ├── CXAUAI_SE1Log.progress.out
│   ├── CXAUAI_SE1ReadsPerGene.out.tab
│   ├── CXAUAI_SE1SJ.out.tab
│   └── CXAUAI_SE1_STARtmp
├── CXAUAI_SE2
│   ├── CXAUAI_SE2Aligned.sortedByCoord.out.bam
│   ├── CXAUAI_SE2Log.final.out
│   ├── CXAUAI_SE2Log.out
│   ├── CXAUAI_SE2Log.progress.out
│   ├── CXAUAI_SE2ReadsPerGene.out.tab
│   ├── CXAUAI_SE2SJ.out.tab
│   └── CXAUAI_SE2_STARtmp
├── CXAUAI_SE3
│   ├── CXAUAI_SE3Aligned.sortedByCoord.out.bam
│   ├── CXAUAI_SE3Log.final.out
│   ├── CXAUAI_SE3Log.out
│   ├── CXAUAI_SE3Log.progress.out
│   ├── CXAUAI_SE3ReadsPerGene.out.tab
│   ├── CXAUAI_SE3SJ.out.tab
│   └── CXAUAI_SE3_STARtmp
├── CXAUAI_SG1
│   ├── CXAUAI_SG1Aligned.sortedByCoord.out.bam
│   ├── CXAUAI_SG1Log.final.out
│   ├── CXAUAI_SG1Log.out
│   ├── CXAUAI_SG1Log.progress.out
│   ├── CXAUAI_SG1ReadsPerGene.out.tab
│   ├── CXAUAI_SG1SJ.out.tab
│   └── CXAUAI_SG1_STARtmp
├── CXAUAI_SG2
│   ├── CXAUAI_SG2Aligned.sortedByCoord.out.bam
│   ├── CXAUAI_SG2Log.final.out
│   ├── CXAUAI_SG2Log.out
│   ├── CXAUAI_SG2Log.progress.out
│   ├── CXAUAI_SG2ReadsPerGene.out.tab
│   ├── CXAUAI_SG2SJ.out.tab
│   └── CXAUAI_SG2_STARtmp
└── CXAUAI_SG3
    ├── CXAUAI_SG3Aligned.sortedByCoord.out.bam
    ├── CXAUAI_SG3Log.final.out
    ├── CXAUAI_SG3Log.out
    ├── CXAUAI_SG3Log.progress.out
    ├── CXAUAI_SG3ReadsPerGene.out.tab
    ├── CXAUAI_SG3SJ.out.tab
    └── CXAUAI_SG3_STARtmp

12 directories, 37 files

```

## 5.Quantify reads 

Environment :
```

mkdir 06.quantify
cd ~/Y_Lipolytica/06.quantify
mkdir alignment-based & mkdir mapping-based
mkdir Ref_transcriptome

[wirawit.s@omic mapping-based]$ tree -L 2
.
├── building_index
│   ├── build_index.sh
│   ├── decoys.txt
│   ├── decoys.txt.bak
│   ├── gentrome.fna
│   ├── parallel.err
│   ├── parallel.log
│   └── salmon_index
├── CXAUAI
│   ├── CXAUAI_n_trim
│   └── CXAUAI_w_trim
└── ESIG
    ├── ESIG_n_trim
    └── ESIG_w_trim
```

8 directories, 6 files


### Mapping-based (Salmon)
##### Create a transcriptome index
I followed [this](https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/) tutorial 

```
cd Ref_transcriptome/
wget \
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/525/GCF_000002525.2_ASM252v1/GCF_000002525.2_ASM252v1_rna.fna.gz
wget \
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/525/GCF_000002525.2_ASM252v1/GCF_000002525.2_ASM252v1_genomic.gtf.gz

```
decompress transcriptome fasta and gtf
change directory to `~/Y_Lipolytica/06.quantify/mapping-based`
make directory `building_index` & change directory to `building_index`

run this script to create decoys and index folder
```sh

#!/bin/bash

module load salmon/0.14.1-mlv2kd4


grep "^>" ~/Y_Lipolytica/05.mapping_reads/Ref_genome/GCF_000002525.2_ASM252v1_genomic.fna \
| cut -d " " -f 1 > decoys.txt

sed -i.bak -e 's/>//g' decoys.txt

cat ~/Y_Lipolytica/06.quantify/Ref_transcriptome/GCF_000002525.2_ASM252v1_rna.fna \
~/Y_Lipolytica/05.mapping_reads/Ref_genome/GCF_000002525.2_ASM252v1_genomic.fna > gentrome.fna

salmon index -t gentrome.fna -i salmon_index --decoys decoys.txt
###############################################

```
result
```
[wirawit.s@omic building_index]$ tree
.
├── build_index.sh
├── decoys.txt
├── decoys.txt.bak
├── gentrome.fna
├── parallel.err
├── parallel.log
└── salmon_index
    ├── duplicate_clusters.tsv
    ├── hash.bin
    ├── header.json
    ├── indexing.log
    ├── quasi_index.log
    ├── refInfo.json
    ├── rsd.bin
    ├── sa.bin
    ├── txpInfo.bin
    └── versionInfo.json

1 directory, 16 files

```

#### Quantifying in mapping-based mode
I created directories (CXAUAI & ESIG )for Salmon outputs in 
```
[wirawit.s@omic mapping-based]$ ls
building_index  CXAUAI  ESIG
[wirawit.s@omic mapping-based]$ tree -L 2
.
├── building_index
│   ├── build_index.sh
│   ├── decoys.txt
│   ├── decoys.txt.bak
│   ├── gentrome.fna
│   ├── parallel.err
│   ├── parallel.log
│   └── salmon_index
├── CXAUAI
│   ├── CXAUAI_n_trim
│   └── CXAUAI_w_trim
└── ESIG
    ├── ESIG_n_trim
    └── ESIG_w_trim


```

example script

```sh

#!/bin/bash

# Set the directory where your data is stored
data_dir=~/Y_Lipolytica/03.trimming/CXAUAI/CXAUAI_S*

module load salmon/0.14.1-mlv2kd4
# Loop through each directory in the data directory
for dir in $data_dir; do
    echo "Salmon is working on files in directory: $dir"
    echo "directory name: ${dir##*/}"
    
    #create directory for output
    mkdir -p ~/Y_Lipolytica/06.quantify/mapping-based/CXAUAI/CXAUAI_n_trim/${dir##*/}
    # Set the file patterns
    pattern_1="$dir/CXAUAI_S*_1.fq.gz"
    pattern_2="$dir/CXAUAI_S*_2.fq.gz"
    
    # List files matching the patterns
    echo "Files matching pattern 1:"
    echo $pattern_1
    
    echo "Files matching pattern 2:"
    echo $pattern_2


    salmon quant \
    -i ~/Y_Lipolytica/06.quantify/mapping-based/building_index/salmon_index \
    -l IU \
    -1 $pattern_1 \
    -2 $pattern_2 \
    --validateMappings \
    -o ~/Y_Lipolytica/06.quantify/mapping-based/CXAUAI/CXAUAI_n_trim/${dir##*/}
    
done

pwd

###############################################
```

#### Result
```
[wirawit.s@omic CXAUAI_n_trim]$ tree -L 2
.
├── CXAUAI_n_trim.sh
├── CXAUAI_SE1
│   ├── aux_info
│   ├── cmd_info.json
│   ├── lib_format_counts.json
│   ├── libParams
│   ├── logs
│   └── quant.sf
├── CXAUAI_SE2
│   ├── aux_info
│   ├── cmd_info.json
│   ├── lib_format_counts.json
│   ├── libParams
│   ├── logs
│   └── quant.sf
├── CXAUAI_SE3
│   ├── aux_info
│   ├── cmd_info.json
│   ├── lib_format_counts.json
│   ├── libParams
│   ├── logs
│   └── quant.sf
├── CXAUAI_SG1
│   ├── aux_info
│   ├── cmd_info.json
│   ├── lib_format_counts.json
│   ├── libParams
│   ├── logs
│   └── quant.sf
├── CXAUAI_SG2
│   ├── aux_info
│   ├── cmd_info.json
│   ├── lib_format_counts.json
│   ├── libParams
│   ├── logs
│   └── quant.sf
├── CXAUAI_SG3
│   ├── aux_info
│   ├── cmd_info.json
│   ├── lib_format_counts.json
│   ├── libParams
│   ├── logs
│   └── quant.sf
├── parallel.err
└── parallel.log

24 directories, 21 files

```

### Alignment-based(Stringtie)

create a directory at 
`~/Y_Lipolytica/06.quantify/count_by_Stringtie/alternative_workflow`

#### Alternative workflow of Stringtie

Stringtie workflow i use is a alternative workflow that does not require "merge" step
because The traditional workflow shown too many MSTRG ids in .tab files
this is the script

```sh
#!/bin/bash

# Iterate over all BAM files
for bam in ~/Y_Lipolytica/05.mapping_reads/CXAUAI/CXAUAI_n_trim/*/*Aligned.sortedByCoord.out.bam; do
    echo "$bam"

    # Extract the prefix using sed and assign it to the variable
    #prefix=$(echo $bam | sed 's|.*/\(CXAUAI_[A-Z0-9]*\)/.*|\1|')
    prefix=$(echo $bam | sed 's|.*/\(CXAUAI_[^/]*\)/.*|\1|')
    echo $prefix

    mkdir ~/Y_Lipolytica/06.quantify/count_by_Stringtie/alternative_workflow/${prefix}
   

    # Run stringtie with the correct output syntax
    ~/Y_Lipolytica/06.quantify/count_by_Stringtie/stringtie-2.2.3/stringtie $bam -p 24 -e -B \
    -o ~/Y_Lipolytica/06.quantify/count_by_Stringtie/alternative_workflow/${prefix}/${prefix}.gtf  \
    -G  ~/Y_Lipolytica/06.quantify/count_by_Stringtie/sol1.gtf \
    -A ~/Y_Lipolytica/06.quantify/count_by_Stringtie/alternative_workflow/${prefix}/${prefix}.tab


done

```

the `sol1.gtf` file it's originate from the script at 

`~/Y_Lipolytica/05.mapping_reads/Ref_genome/GCF_000002525.2_ASM252v1_genomic.gtf`

with a awk command because when i use 
`GCF_000002525.2_ASM252v1_genomic.gtf`
it' shown Error: no valid ID found for GFF record

```sh
awk '$3 != "gene" ' my_annotation.gtf > my_annotation_no_genes.gtf
```

#### Result 

```
[wirawit.s@master alternative_workflow]$ tree -L 1
.
├── CXAUAI_SE1
├── CXAUAI_SE1_test.gtf
├── CXAUAI_SE1_test.tab
├── CXAUAI_SE2
├── CXAUAI_SE3
├── CXAUAI_SG1
├── CXAUAI_SG2
├── CXAUAI_SG2_test.gtf
├── CXAUAI_SG2_test.tab
├── CXAUAI_SG3
├── ESIG_SE-1
├── ESIG_SE-2
├── ESIG_SE-3
├── ESIG_SG-1
├── ESIG_SG-2
├── ESIG_SG-3
├── Stringtie_CXAUAI.sh
├── Stringtie_ESIG.sh
├── Stringtie_test.sh
├── test
├── tmp2
└── TPM_all

14 directories, 8 files
```

for example
```

[wirawit.s@master alternative_workflow]$ ls CXAUAI_SE1
CXAUAI_SE1.gtf  CXAUAI_SE1.tab  e2t.ctab  e_data.ctab  i2t.ctab  i_data.ctab  t_data.ctab

```

## 6.Create count matrix(s)
### Create a raw count matrix from HTSeq 
When you use STAR and specified --quant mode it's generate  `ReadPerGene.out.tab`

```
[wirawit.s@master CXAUAI_SE1]$ ls
CXAUAI_SE1Aligned.sortedByCoord.out.bam  CXAUAI_SE1Log.out           CXAUAI_SE1ReadsPerGene.out.tab  CXAUAI_SE1_STARtmp
CXAUAI_SE1Log.final.out                  CXAUAI_SE1Log.progress.out  CXAUAI_SE1SJ.out.tab

#the inside look like this
[wirawit.s@master CXAUAI_SE1]$ head CXAUAI_SE1ReadsPerGene.out.tab
N_unmapped      3308648 3308648 3308648
N_multimapping  961183  961183  961183
N_noFeature     4959885 39271256        40404876
N_ambiguous     66147   17284   17223
YALI0_A00110g   244     115     129
YALI0_A00132g   134247  68093   66174
YALI0_A00154g   22980   11667   11315
YALI0_A00176g   2289    1085    1204
YALI0_A00198g   187     88      99
YALI0_A00212g   822     415     407

```
which explain in the HTSeq 
"For `stranded=no` (2nd column), a read is considered overlapping with a feature regardless of whether it is mapped to the same or the opposite strand as the feature. For `stranded=yes` (3rd column) and single-end reads, the read has to be mapped to the same strand as the feature. For paired-end reads, the first read has to be on the same strand and the second read on the opposite strand. For `stranded=reverse` (4th column), these rules are reversed."

i copied all `ReadPerGene.out.tab` to a folder then i use a script to create count matrix 

```sh
#!/bin/bash

paste *ReadsPerGene.out.tab | grep -v "N" |awk '{printf "%s\t", $1}{for (i=2;i<=NF;i+=4) printf "%s\t", $i; printf "\n" }' > tmp2


sed -e "1igene_name\t$(ls *ReadsPerGene.out.tab | sed 's/ReadsPerGene.out.tab//g' | tr '\n' '\t')" tmp2 > star_raw_counts_matrix_v2.txt

```
#### Result
result in a tab delimited raw count file which a be use for further analysis
```
[wirawit.s@master ReadPerGene_sum]$ head star_raw_counts_matrix_v2.txt
gene_name       CXAUAI_SE1      CXAUAI_SE2      CXAUAI_SE3      CXAUAI_SG1      CXAUAI_SG2      CXAUAI_SG3      ESIG_SE-1       ESIG_SE-2       ESIG_SE-3       ESIG_SG-1       ESIG_SG-2       ESIG_SG-3       trim_CXAUAI_SE1 trim_CXAUAI_SE2   trim_CXAUAI_SE3 trim_CXAUAI_SG1 trim_CXAUAI_SG2 trim_CXAUAI_SG3 trim_ESIG_SE-1  trim_ESIG_SE-2  trim_ESIG_SE-3  trim_ESIG_SG-1  trim_ESIG_SG-2  trim_ESIG_SG-3
YALI0_A00110g   244     201     142     109     83      75      1849    79      91      280     39      36      231     191     134     107     83      74      1845    79      90      275     38      36
YALI0_A00132g   134247  100109  110260  107963  103803  116114  61260   74398   50562   91359   52590   65402   130714  96646   107215  105697  102718  113990  61058   73970   50284   88978   52434   64972

```
### Create a raw count matrix from Salmon `quant.sf` 
Similar to STAR  `ReadPerGene.out.tab` i copied and have `quant.sf`  in one folder
```
[wirawit.s@master quant_sf_sum]$ ls
CXAUAI_SE1_quant.sf  CXAUAI_SG2_quant.sf  ESIG_SE-3_quant.sf  raw_sal_counts_matrix.txt  trim_CXAUAI_SE1_quant.sf  trim_CXAUAI_SG2_quant.sf  trim_ESIG_SE-3_quant.sf
CXAUAI_SE2_quant.sf  CXAUAI_SG3_quant.sf  ESIG_SG-1_quant.sf  script.sh                  trim_CXAUAI_SE2_quant.sf  trim_CXAUAI_SG3_quant.sf  trim_ESIG_SG-1_quant.sf
CXAUAI_SE3_quant.sf  ESIG_SE-1_quant.sf   ESIG_SG-2_quant.sf  Test                       trim_CXAUAI_SE3_quant.sf  trim_ESIG_SE-1_quant.sf   trim_ESIG_SG-2_quant.sf
CXAUAI_SG1_quant.sf  ESIG_SE-2_quant.sf   ESIG_SG-3_quant.sf  tmp                        trim_CXAUAI_SG1_quant.sf  trim_ESIG_SE-2_quant.sf   trim_ESIG_SG-3_quant.sf

#for example 
[wirawit.s@master quant_sf_sum]$ head CXAUAI_SE1_quant.sf
Name    Length  EffectiveLength TPM     NumReads
XM_002142944.1  309     67.461  93.529084       630.000
XM_002142945.1  1893    1639.348        3.964917        649.000
XM_002142946.1  168     9.480   123.602407      117.000
XM_002142947.1  165     8.328   22.848395       19.000
XM_002142948.1  1248    994.348 0.090649        9.000
XM_002142949.1  288     53.596  24.479328       131.000
XM_002142950.1  534     280.396 329.928679      9237.000
XM_002142951.1  78      79.000  0.000000        0.000
XM_002142952.1  3294    3040.348        28.134944       8541.000

```
and run this script
```sh
[wirawit.s@master quant_sf_sum]$ cat script.sh
#!/bin/bash

paste *quant.sf | grep -v "N" |awk '{printf "%s\t", $1}{for (i=5;i<=NF;i+=5) printf "%s\t", $i; printf "\n" }' > tmp


sed -e "1itranscript_name\t$(ls *quant.sf | sed 's/quant.sf//g' | tr '\n' '\t')" tmp > raw_sal_counts_matrix.txt

```

#### Result
```
[wirawit.s@master quant_sf_sum]$ head raw_sal_counts_matrix.txt
transcript_name CXAUAI_SE1_     CXAUAI_SE2_     CXAUAI_SE3_     CXAUAI_SG1_     CXAUAI_SG2_     CXAUAI_SG3_     ESIG_SE-1_      ESIG_SE-2_      ESIG_SE-3_      ESIG_SG-1_      ESIG_SG-2_      ESIG_SG-3_      trim_CXAUAI_SE1_        trim_CXAUAI_SE2_  trim_CXAUAI_SE3_        trim_CXAUAI_SG1_        trim_CXAUAI_SG2_        trim_CXAUAI_SG3_        trim_ESIG_SE-1_ trim_ESIG_SE-2_ trim_ESIG_SE-3_ trim_ESIG_SG-1_ trim_ESIG_SG-2_ trim_ESIG_SG-3_
XM_002142944.1  630.000 301.000 373.000 562.000 411.000 409.000 156.000 111.000 225.000 363.000 274.000 297.000 874.000 459.000 497.000 805.000 610.000 592.000 236.000 187.000 327.000 527.000 425.000 499.000
XM_002142945.1  649.000 472.000 622.000 628.000 504.000 586.000 237.000 242.000 257.000 225.000 106.000 154.000 658.000 472.000 627.000 636.000 518.000 593.000 248.000 250.000 263.000 226.000 106.000 157.000

```

### Create a TPM matrix from Stringtie `.tab` files 
```
[wirawit.s@master TPM_all]$ ls
CXAUAI_SE1.tab  CXAUAI_SE3.tab  CXAUAI_SG2.tab  ESIG_SE-1.tab  ESIG_SE-3.tab  ESIG_SG-2.tab  script.sh  TPM_sum.txt
CXAUAI_SE2.tab  CXAUAI_SG1.tab  CXAUAI_SG3.tab  ESIG_SE-2.tab  ESIG_SG-1.tab  ESIG_SG-3.tab  tmp

#example

[wirawit.s@master TPM_all]$ head CXAUAI_SE1.tab
Gene ID Gene Name       Reference       Strand  Start   End     Coverage        FPKM    TPM
YALI0_A00110g   -       NC_006067.1     +       2659    5277    22.220354       1.071615        1.274624
YALI0_A00132g   -       NC_006067.1     +       7045    8880    18272.744141    881.234375      1048.177612
YALI0_A00154g   -       NC_006067.1     +       11559   12653   5057.270996     243.895569      290.099731

```
script 

```sh
#!/bin/bash

paste *.tab | grep -v "Gene" |awk '{printf "%s\t", $1}{for (i=9;i<=NF;i+=9) printf "%s\t", $i; printf "\n" }' > tmp


sed -e "1igene_name\t$(ls *.tab | sed 's/.tab//g' | tr '\n' '\t')" tmp > TPM_sum.txt

```
#### Result 
```
[wirawit.s@master TPM_all]$ head TPM_sum.txt
gene_name       CXAUAI_SE1      CXAUAI_SE2      CXAUAI_SE3      CXAUAI_SG1      CXAUAI_SG2      CXAUAI_SG3      ESIG_SE-1       ESIG_SE-2       ESIG_SE-3   ESIG_SG-1       ESIG_SG-2       ESIG_SG-3
YALI0_A00110g   1.274624        1.474144        0.911001        0.723814        0.591808        0.481355        22.867292       0.938317        1.115107    2.669852        0.464589        0.411186

```
### Create a raw count matrix with `prepDE.py`
after obtain Stringtie output
```
[wirawit.s@master outputs]$ tree -L 2
.
├── CXAUAI_SE1
│   ├── CXAUAI_SE1.gtf
│   ├── CXAUAI_SE1.tab
│   ├── e2t.ctab
│   ├── e_data.ctab
│   ├── i2t.ctab
│   ├── i_data.ctab
│   └── t_data.ctab
├── CXAUAI_SE2
│   ├── CXAUAI_SE2.gtf
│   ├── CXAUAI_SE2.tab
│   ├── e2t.ctab
│   ├── e_data.ctab
│   ├── i2t.ctab
│   ├── i_data.ctab
│   └── t_data.ctab
├── CXAUAI_SE3
│   ├── CXAUAI_SE3.gtf
│   ├── CXAUAI_SE3.tab
│   ├── e2t.ctab
│   ├── e_data.ctab
│   ├── i2t.ctab
│   ├── i_data.ctab
│   └── t_data.ctab
├── CXAUAI_SG1
│   ├── CXAUAI_SG1.gtf
│   ├── CXAUAI_SG1.tab
│   ├── e2t.ctab
│   ├── e_data.ctab
│   ├── i2t.ctab
│   ├── i_data.ctab
│   └── t_data.ctab
├── CXAUAI_SG2
│   ├── CXAUAI_SG2.gtf
│   ├── CXAUAI_SG2.tab
│   ├── e2t.ctab
│   ├── e_data.ctab
│   ├── i2t.ctab
│   ├── i_data.ctab
│   └── t_data.ctab
├── CXAUAI_SG3
│   ├── CXAUAI_SG3.gtf
│   ├── CXAUAI_SG3.tab
│   ├── e2t.ctab
│   ├── e_data.ctab
│   ├── i2t.ctab
│   ├── i_data.ctab
│   └── t_data.ctab
├── ESIG_SE-1
│   ├── e2t.ctab
│   ├── e_data.ctab
│   ├── ESIG_SE-1.gtf
│   ├── ESIG_SE-1.tab
│   ├── i2t.ctab
│   ├── i_data.ctab
│   └── t_data.ctab
├── ESIG_SE-2
│   ├── e2t.ctab
│   ├── e_data.ctab
│   ├── ESIG_SE-2.gtf
│   ├── ESIG_SE-2.tab
│   ├── i2t.ctab
│   ├── i_data.ctab
│   └── t_data.ctab
├── ESIG_SE-3
│   ├── e2t.ctab
│   ├── e_data.ctab
│   ├── ESIG_SE-3.gtf
│   ├── ESIG_SE-3.tab
│   ├── i2t.ctab
│   ├── i_data.ctab
│   └── t_data.ctab
├── ESIG_SG-1
│   ├── e2t.ctab
│   ├── e_data.ctab
│   ├── ESIG_SG-1.gtf
│   ├── ESIG_SG-1.tab
│   ├── i2t.ctab
│   ├── i_data.ctab
│   └── t_data.ctab
├── ESIG_SG-2
│   ├── e2t.ctab
│   ├── e_data.ctab
│   ├── ESIG_SG-2.gtf
│   ├── ESIG_SG-2.tab
│   ├── i2t.ctab
│   ├── i_data.ctab
│   └── t_data.ctab
├── ESIG_SG-3
│   ├── e2t.ctab
│   ├── e_data.ctab
│   ├── ESIG_SG-3.gtf
│   ├── ESIG_SG-3.tab
│   ├── i2t.ctab
│   ├── i_data.ctab
│   └── t_data.ctab
```

then i run this script
```sh
[wirawit.s@master outputs]$ cat prep_forDE.sh
#!/bin/bash
~/Y_Lipolytica/06.quantify/count_by_Stringtie/alternative_workflow_v_221/stringtie-2.2.1/prepDE.py3
```
#### Result
these are result from
```
[wirawit.s@master outputs]$ head *csv
==> gene_count_matrix.csv <==
gene_id,CXAUAI_SE1,CXAUAI_SE2,CXAUAI_SE3,CXAUAI_SG1,CXAUAI_SG2,CXAUAI_SG3,ESIG_SE-1,ESIG_SE-2,ESIG_SE-3,ESIG_SG-1,ESIG_SG-2,ESIG_SG-3
YALI0_A00715r,0,0,0,4,0,0,0,0,0,0,0,0
YALI0_A01287r,0,6,3,5,4,8,0,0,1,0,0,0
YALI0_A01914r,0,0,0,0,0,0,0,0,0,0,0,0
YALI0_A03047r,0,0,0,0,0,0,0,0,0,0,1,1
YALI0_A03157r,0,0,0,0,2,1,2,0,0,0,0,0
YALI0_A04389r,0,0,0,0,0,4,0,0,0,0,0,0
YALI0_A04763r,0,4,1,0,0,0,1,0,2,0,2,0
YALI0_A04829r,0,0,0,0,0,0,0,0,0,0,0,0
YALI0_A05049r,0,0,3,0,2,2,0,0,0,0,0,0

==> transcript_count_matrix.csv <==
transcript_id,CXAUAI_SE1,CXAUAI_SE2,CXAUAI_SE3,CXAUAI_SG1,CXAUAI_SG2,CXAUAI_SG3,ESIG_SE-1,ESIG_SE-2,ESIG_SE-3,ESIG_SG-1,ESIG_SG-2,ESIG_SG-3
unassigned_transcript_3,0,0,0,4,0,0,0,0,0,0,0,0
unassigned_transcript_6,0,6,3,5,4,8,0,0,1,0,0,0
unassigned_transcript_9,0,0,0,0,0,0,0,0,0,0,0,0
unassigned_transcript_12,0,0,0,0,0,0,0,0,0,0,1,1
unassigned_transcript_13,0,0,0,0,2,1,2,0,0,0,0,0
unassigned_transcript_24,0,0,0,0,0,4,0,0,0,0,0,0
unassigned_transcript_29,0,4,1,0,0,0,1,0,2,0,2,0
unassigned_transcript_30,0,0,0,0,0,0,0,0,0,0,0,0
unassigned_transcript_31,0,0,3,0,2,2,0,0,0,0,0,0

```

# And for R programing
## Plot from TPM
Boxplot and PCA
```r
setwd("TPM/")
TPM.tab <- read.delim("TPM_sum_v2.txt")

TPM.tab$X <- NULL
rownames(TPM.tab) <- TPM.tab$gene_name
TPM.tab$gene_name <- NULL

TPM.matrix<- as.matrix(TPM.tab)

par(mfrow = c(1,1))

boxplot(log2(TPM.matrix),ylab=expression('Log'[2]~'TPM'),las=2,main="TPM counts")


#find gene(row) with no variable (all 1) every sample(col) 
novar <- as.numeric(which(apply(TPM.matrix, 1, var)==0))

TPM.matrix.allvar <- TPM.matrix[-novar,]

#double chech for non variable genes
(which(apply(TPM.matrix.allvar, 1, var)==0))

#this will not work
TPM.pca <- prcomp(t(TPM.matrix), scale = TRUE) 

TPM.pca <- prcomp(t(TPM.matrix.allvar), scale = TRUE) 

# summary of the  
# prcomp object 
summary(TPM.pca)

plot(TPM.pca$x[,1], TPM.pca$x[,2])


pca.var <- TPM.pca$sdev^2
pca.var.perc <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.perc, main = "Scree Plot",
        xlab = "Principal Component",
        ylab = "Percent Variation")



#prep data for ggplot
pca.Data <- data.frame(Samle = rownames(TPM.pca$x),
                       X = TPM.pca$x[,1],
                       Y = TPM.pca$x[,2])
library(ggplot2)
#use pca.Data tell that x is X in pca.Data y is Y in pca.Data
#and label is Sample in pca.Data which is row name from TPM.pca$x
ggplot(data = pca.Data,aes(x = X , y = Y, label=Samle)) +
         geom_text() +
  #name x and y label using pca.var.perc
         xlab(paste("PC1 - ",pca.var.perc[1], "%", sep = "")) +
         ylab(paste("PC2 - ",pca.var.perc[2], "%", sep = "")) +
         theme_light()+
         ggtitle("PCA Graph")


```

TPM Heatmap 
```r
TPM <- read.csv("TPM_sum_v2.txt",
                sep = "\t",
                header = TRUE)
TPM$X <- NULL

rownames(TPM) <- TPM$gene_name

TPM$gene_name <- NULL



TPM.matrix <- as.matrix(TPM)

TPM.matrixplus.one <- TPM.matrix + 1

matrixlog.t <- log(TPM.matrixplus.one, base = 2)


z_new <- t(scale(t(TPM.matrix)))

heatmap(TPM.matrix,
        scale = "row")

z_new <- na.omit(z_new)

heatmap(z_new,
        scale = "row",
        col = my_colors(100))

my_colors<- colorRampPalette(c("green", "white", "red"))

install.packages('pheatmap') 
library(pheatmap)
pheatmap(z_new,
         fontsize = 9,
         cutree_rows = 5,
         cutree_cols = 2)

```

DESeq2
