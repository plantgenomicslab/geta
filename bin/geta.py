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

genome = os.path.abspath(genome)
if not genome:
    raise ValueError("No genome fasta input")

protein = os.path.abspath(protein)
if not ((pe1 and pe2) or single_end or protein):
    raise ValueError("No RNA-Seq short reads or homologous proteins as input")

if not (augustus_species or use_existed_augustus_species):
    raise ValueError("No Augustus species provided")

if use_existed_augustus_species:
    species_config_dir = os.environ["AUGUSTUS_CONFIG_PATH"]
    species_config_dir = os.path.join(species_config_dir, "species", use_existed_augustus_species)
    if not os.path.exists(species_config_dir):
        raise ValueError("The AUGUSUTS HMM files of {} does not exists!".format(use_existed_augustus_species))
    else:
        augustus_species = use_existed_augustus_species

if RM_lib:
    RM_lib = os.path.abspath(RM_lib)

if config:
    config = os.path.abspath(config)

pe_reads = {}
se_reads = {}

if pe1 and pe2:
    pe1 = pe1.split(",")
    pe2 = pe2.split(",")
    if len(pe1) != len(pe2):
        raise ValueError("the input file number of -1 was not equal to -2.")
    for i in range(len(pe1)):
        pe1[i] = os.path.abspath(pe1[i])
        pe2[i] = os.path.abspath(pe2[i])
        pe_file = pe1[i] + "\t" + pe2[i]
        pe_reads[pe_file] = 1

if single_end:
    se = single_end.split(",")
    for i in range(len(se)):
        se[i] = os.path.abspath(se[i])
        se_reads[se[i]] = 1
   
out_prefix = out_prefix or "out"
cpu = cpu or 4

HMM_db = {}
if HMM_db:
    for db in HMM_db.split(","):
        db = os.path.abspath(db)
        HMM_db[db] = os.path.basename(db)

BLASTP_db = {}
if BLASTP_db:
    for db in BLASTP_db.split(","):
        BLASTP_db[db] = os.path.basename(db)


config = {
    'RepeatMasker': '-e ncbi -gff',
    'trimmomatic': 'TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 TOPHRED33',
    'hisat2-build': '-p 1',
    'hisat2': '--min-intronlen 20 --max-intronlen 20000 --dta --score-min L,0.0,-0.4',
    'sam2transfrag': '--fraction 0.05 --min_expressed_base_depth 2 --max_expressed_base_depth 50 --min_junction_depth 2 --max_junction_depth 50 --min_fragment_count_per_transfrags 10 --min_intron_length 20',
    'TransDecoder.LongOrfs': '-m 100 -G universal',
    'TransDecoder.Predict': '--retain_long_orfs_mode dynamic',
    'homolog_genewise': '--coverage_ratio 0.4 --evalue 1e-9',
    'homolog_genewiseGFF2GFF3': '--min_score 15 --gene_prefix genewise --filterMiddleStopCodon',
    'geneModels2AugusutsTrainingInput': '--min_evalue 1e-9 --min_identity 0.8 --min_coverage_ratio 0.8 --min_cds_num 2 --min_cds_length 450 --min_cds_exon_ratio 0.60',
    'BGM2AT': '--min_gene_number_for_augustus_training 500 --gene_number_for_accuracy_detection 200 --min_gene_number_of_optimize_augustus_chunk 50 --max_gene_number_of_optimize_augustus_chunk 200',
    'prepareAugusutusHints': '--margin 20',
    'paraAugusutusWithHints': '--gene_prefix augustus --min_intron_len 20',
    'paraCombineGeneModels': '--overlap 30 --min_augustus_transcriptSupport_percentage 10.0 --min_augustus_intronSupport_number 1 --min_augustus_intronSupport_ratio 0.01',
    'pickout_better_geneModels_from_evidence': '--overlap_ratio 0.2 --ratio1 2 --ratio2 1.5 --ratio3 0.85 --ratio4 0.85',
    'fillingEndsOfGeneModels': '--start_codon ATG --stop_codon TAG,TGA,TAA',
    'alternative_splicing_analysis': '--min_intron_depth 1 --min_base_depth_ratio_for_ref_specific_intron 0.3 --min_intron_depth_ratio_for_evidence_specific_intron 0.2 --min_base_depth_ratio_for_common_intron 0.2 --min_gene_depth 10 --min_transcript_confidence_for_output 0.05 --transcript_num_for_output_when_all_low_confidence 8 --added_mRNA_ID_prefix t',
    'GFF3_extract_TranscriptID_for_filtering': '--min_CDS_ratio 0.3 --min_CDS_length 600 --max_repeat_overlap_ratio 0.3 --ignore_repeat_Name Simple_repeat,Low_complexity,Satellite,Unknown,Tandem_repeat',
    'para_hmmscan': '--evalue1 1e-5 --evalue2 1e-3 --hmm_length 80 --coverage 0.25 --no_cut_ga --chunk 20 --hmmscan_cpu 2',
    'diamond': '--sensitive --max-target-seqs 20 --evalue 1e-5 --id 10 --index-chunks 1 --block-size 5',
    'parsing_blast_result.pl': '--evalue 1e-9 --identity 0.8 --coverage 0.8 --min_length 200 --min_cds_num 2 --min_cds_length 450 --min_cds_exon_ratio 0.60'
}



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