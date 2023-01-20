#!/usr/bin/env python
"""
Geta.py
"""

import os
import argparse
import subprocess
import datetime
import shutil

parser = argparse.ArgumentParser()
args = parser.parse_args()

script_dir = os.path.dirname(os.path.abspath(__file__))
script_dir = script_dir.rstrip("/bin")

parser = argparse.ArgumentParser(description='description of your script')
parser.add_argument('--RM_species', type=str, required=True, help='description of RM_species')
parser.add_argument('--RM_lib', type=str, help='description of RM_lib')
parser.add_argument('--genome', type=str, required=True, help='description of genome')
parser.add_argument('--out_prefix', type=str, required=True, help='description of out_prefix')
parser.add_argument('-1', type=str, help='description of pe1')
parser.add_argument('-2', type=str, help='description of pe2')
parser.add_argument('-S', type=str, help='description of single_end')
parser.add_argument('--protein', type=str, help='description of protein')
parser.add_argument('--cpu', type=int, help='description of cpu')
parser.add_argument('--strand_specific', action='store_true', help='description of strand_specific')
parser.add_argument('--augustus_species', type=str, help='description of augustus_species')
parser.add_argument('--use_existed_augustus_species', type=str, help='description of use_existed_augustus_species')
parser.add_argument('--HMM_db', type=str, help='description of HMM_db')
parser.add_argument('--BLASTP_db', type=str, help='description of BLASTP_db')
parser.add_argument('--gene_prefix', type=str, help='description of gene_prefix')
parser.add_argument('--enable_augustus_training_iteration', action='store_true', help='description of enable_augustus_training_iteration')
parser.add_argument('--config', type=str, help='description of config')
parser.add_argument('--augustus_species_start_from', type=str, help='description of augustus_species_start_from')
parser.add_argument('--no_alternative_splicing_analysis', action='store_true', help='description of no_alternative_splicing_analysis')

args = parser.parse_args()

software_list = ["RepeatMasker","RepeatModeler","samtools","java","blastn","hisat2","hmmscan"]
check_software_existence(software_list)

#check_software_existence("RepeatMasker")
#check_software_existence("samtools")

#def check_software_existence(software_name:str)->None:
def check_software_existence(software_name:list)->None:
    for software in software_list:
        if shutil.which(software_name) is not None:
            print(f"{software_name} is installed")
        else:
            print(f"{software_name} is not installed")



# Get current working directory
pwd = subprocess.check_output(["pwd"]).strip().decode("utf-8")
print("PWD:", pwd)

# Get current date and time
print(datetime.datetime.now(), ": CMD:", " ".join(sys.argv))






"""
Usage:
    perl $0 [options]

For example:
    perl $0 --RM_species Embryophyta --genome genome.fasta -1 liba.1.fq.gz,libb.1.fq.gz -2 liba.2.fq.gz,libb.2.fq.gz --protein homolog.fasta --augustus_species oryza_sativa_20171120 --out_prefix out --config conf.txt --cpu 80 --gene_prefix OS01Gene --pfam_db ~/opt/biosoft/hmmer-3.1b2/Pfam-AB.hmm

Parameters:
[General]
    --genome <string>     Required
    genome file in fasta format.

    -1 <string> -2 <string>    Not Required but Recommened
    fastq format files contain of paired-end RNA-seq data. if you have data come from multi librarys, input multi fastq files separated by comma. the compress file format .gz also can be accepted.

    -S <string>    Not Required, a option when -1 and -2 were not provided
    fastq format file contains of single-end RNA-seq data. if you have data come from multi librarys, input multi fastq files separated by comma. the compress file format .gz also can be accepted.

    --protein <string>    Required
    homologous protein sequences (derived from multiple species would be recommended) file in fasta format.

    --augustus_species <string>    Required when --use_existed_augustus_species were not provided
    species identifier for Augustus. the relative hmm files of augustus training will be created with this prefix. if the relative hmm files of augustus training exists, the program will delete the hmm files directory firstly, and then start the augustus training steps.

    --use_existed_augustus_species <string>    Required when --augustus_species were not provided
    species identifier for Augustus. This parameter is conflict with --augustus_species. When this parameter set, the --augustus_species parameter will be invalid, and the relative hmm files of augustus training should exists, and the augustus training step will be skipped (this will save lots of runing time).

[other]
    --out_prefix <string>    default: out
    the prefix of outputs.

    --config <string>    default: None
    Input a file containing the parameters of several main programs (such as trimmomatic, hisat2 and augustus) during the pipeline. If you do not input this file, the default parameters should be suitable for most situation.
    
    --RM_species <string>    default: None
    species identifier for RepeatMasker. The acceptable value of this parameter can be found in file $dirname/RepeatMasker_species.txt. Such as, Eukaryota for eucaryon, Fungi for fungi, Viridiplantae for plants, Metazoa for animals. The repeats in genome sequences would be searched aganist the Repbase database when this parameter set. 

    --RM_lib <string>    default: None
    A fasta file of repeat sequences. Generally to be the result of RepeatModeler. If not set, RepeatModeler will be used to product this file automaticly, which shall time-consuming.

    --augustus_species_start_from <string>    default: None
    species identifier for Augustus. The optimization step of Augustus training will start from the parameter file of this species, so it may save much time when setting a close species.

    --cpu <int>    default: 4
    the number of threads.

    --strand_specific    default: False
    enable the ability of analysing the strand-specific information provided by the tag "XS" from SAM format alignments. If this parameter was set, the paramter "--rna-strandness" of hisat2 should be set to "RF" usually.

    --HMM_db <string>    default: None
    the absolute path of protein family HMM database which was used for filtering of false positive gene models. multiple databases can be input, and the prefix of database files should be seperated by comma.

    --BLASTP_db <string>    default: None
    the absolute path of protein family diamond database which was used for filtering of false positive gene models. 若该参数没有设置，程序会以homologous protein构建diamond数据库，进行基因模型过滤。multiple databases can be input, and the prefix of database files should be seperated by comma.

    --gene_prefix <string>    default: gene
    the prefix of gene id shown in output file.

    --enable_augustus_training_iteration    default: False
    开启augustus_training_iteration，运行在第一次Augustus training后，根据基因预测的结果，选择有证据支持的基因模型，再一次进行Augustus training（迭代）。此举会消耗较多计算时间，且可能对基因预测没有改进，或产生不好的影响。

    --no_alternative_splicing_analysis    default: None
    添加该参数后，程序不会进行可变剪接分析。


This script was tested on CentOS 8.4 with such softwares can be run directly in terminal:
01. ParaFly
02. java (version: 1.8.0_282)
03. hisat2 (version: 2.1.0)
04. samtools (version: 1.10)
05. hmmscan (version: 3.3.1)
06. makeblastdb/tblastn/blastp (version: 2.6.0)
07. RepeatMasker (version: 4.1.2-p1)
08. RepeatModeler (version: 2.0.3)
09. genewise (version: 2.4.1)
10. augustus/etraining (version: 3.4.0)
11. diamond (version 2.0.2.140)

Version: 2.5.6

USAGE
"""