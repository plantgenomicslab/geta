rule combine_gene_models:
    input:
        augustus = "5.augustus/augustus.gff3",
        transcript = "3.transcript/transfrag.genome.gff3",
        homolog = "4.homolog/genewise.gff3",
        hints = "5.augustus/hints.gff3",
        config = "{params.dirname}/bin/paraCombineGeneModels",
        genome = "{params.genome}"
    output:
        protected("6.combineGeneModels/01.paraCombineGeneModels.ok")
    params:
        dirname = "/path/to/dirname", # Replace with the correct path
        cpu = 8, # Replace with desired number of CPUs
        genome = "/path/to/genome" # Replace with the correct path
    shell:
        """
        mkdir -p 6.combineGeneModels
        cd 6.combineGeneModels
        pwd > /dev/null
        echo "Creating geneModels.Readme..."
        echo '...content...' > geneModels.Readme # Replace '...content...' with the original Readme content
        echo "Performing first round of gene prediction result integration..."
        {input.config} --cpu {params.cpu} {input.augustus} {input.transcript} {input.homolog} {input.hints} > /dev/null
        echo "Running GFF3Clear on combine.1.gff3 and combine.2.gff3..."
        {params.dirname}/bin/GFF3Clear --genome {input.genome} --no_attr_add --coverage 0.8 combine.1.gff3 > geneModels.a.gff3 2> /dev/null
        {params.dirname}/bin/GFF3Clear --genome {input.genome} --no_attr_add --coverage 0.8 combine.2.gff3 > geneModels.b.gff3 2> /dev/null
        echo "Performing Perl replacement on geneModels.a.gff3 and geneModels.b.gff3..."
        perl -p -i -e 's/(=[^;]+)\\.t1/\\$1.t01/g' geneModels.a.gff3 geneModels.b.gff3
        echo "Touching 01.paraCombineGeneModels.ok..."
        touch 01.paraCombineGeneModels.ok
        """
rule second_round_gene_prediction:
    input:
        gff3 = "5.augustus/training/geneModels.gff3",
        geneModels_a = "6.combineGeneModels/geneModels.a.gff3",
        config = "{params.dirname}/bin/pickout_better_geneModels_from_evidence",
        genome = "{params.genome}"
    output:
        protected("6.combineGeneModels/02.pickout_better_geneModels_from_evidence.ok")
    params:
        dirname = "/path/to/dirname" # Replace with the correct path
    shell:
        """
        cd 6.combineGeneModels
        echo "Performing second round of gene prediction result integration..."
        perl -p -e 's/(=[^;]+)\\.t1/\\$1.t01/g;' {input.gff3} > geneModels.c.gff3
        {input.config} geneModels.a.gff3 geneModels.c.gff3 > picked_evidence_geneModels.gff3 2> picked_evidence_geneModels.log
        {params.dirname}/bin/GFF3Clear --genome {input.genome} --no_attr_add picked_evidence_geneModels.gff3 geneModels.a.gff3 > geneModels.d.gff3 2> /dev/null
        perl -p -i -e 's/Integrity=[^;]+;?//g' geneModels.d.gff3
        touch 02.pickout_better_geneModels_from_evidence.ok
        """

rule filling_ends_gene_models:
    input:
        config = "{params.dirname}/bin/fillingEndsOfGeneModels",
        genome = "{params.genome}",
        geneModels_f = "6.combineGeneModels/geneModels.f.gff3",
        geneModels_d = "6.combineGeneModels/geneModels.d.gff3"
    output:
        protected("6.combineGeneModels/03.fillingEndsOfGeneModels.ok")
    params:
        dirname = "/path/to/dirname" # Replace with the correct path
    shell:
        """
        cd 6.combineGeneModels
        echo "Filling ends of gene models..."
        {input.config} --filling_need_transcriptID filling_need_transcriptID.txt --nonCompletedGeneModels {input.geneModels_f} {input.genome} {input.geneModels_d} > geneModels.e.gff3 2> fillingEndsOfGeneModels.1.log
        touch 03.fillingEndsOfGeneModels.ok
        """
rule alternative_splicing_analysis:
    input:
        config = "{params.dirname}/bin/paraAlternative_splicing_analysis",
        genome = "{params.genome}",
        geneModels_b = "6.combineGeneModels/geneModels.b.gff3",
        geneModels_e = "6.combineGeneModels/geneModels.e.gff3",
        geneModels_f = "6.combineGeneModels/geneModels.f.gff3",
        intron = "3.transcript/intron.txt",
        base_depth = "3.transcript/base_depth.txt"
    output:
        protected("6.combineGeneModels/04.alternative_splicing_analysis.ok")
    params:
        dirname = "/path/to/dirname", # Replace with the correct path
        cpu = 8 # Adjust according to your system
    shell:
        """
        cd 6.combineGeneModels
        if [[ ! -v no_alternative_splicing_analysis ]]; then
            echo "Performing alternative splicing analysis..."
            {input.config} --tmp_dir paraAlternative_splicing_analysis.gb.tmp --cpu {params.cpu} {input.geneModels_b} {input.intron} {input.base_depth} > geneModels.gb_AS.gff3 2> alternative_splicing.gb.stats && {params.dirname}/bin/GFF3_add_CDS_for_transcript {input.genome} geneModels.gb_AS.gff3 > geneModels.gb.gff3
            {input.config} --tmp_dir paraAlternative_splicing_analysis.ge.tmp --cpu {params.cpu} {input.geneModels_e} {input.intron} {input.base_depth} > geneModels.ge_AS.gff3 2> alternative_splicing.ge.stats && {params.dirname}/bin/GFF3_add_CDS_for_transcript {input.genome} geneModels.ge_AS.gff3 > geneModels.ge.gff3
            {input.config} --tmp_dir paraAlternative_splicing_analysis.gf.tmp --cpu {params.cpu} {input.geneModels_f} {input.intron} {input.base_depth} > geneModels.gf_AS.gff3 2> alternative_splicing.gf.stats && {params.dirname}/bin/GFF3_add_CDS_for_transcript {input.genome} geneModels.gf_AS.gff3 > geneModels.gf.gff3
        else
            echo "Adding CDS for transcripts without alternative splicing analysis..."
            {params.dirname}/bin/GFF3_add_CDS_for_transcript {input.genome} {input.geneModels_b} > geneModels.gb.gff3
            {params.dirname}/bin/GFF3_add_CDS_for_transcript {input.genome} {input.geneModels_e} > geneModels.ge.gff3
            {params.dirname}/bin/GFF3_add_CDS_for_transcript {input.genome} {input.geneModels_f} > geneModels.gf.gff3
        fi
        touch 04.alternative_splicing_analysis.ok
        """
rule extract_TranscriptID_for_filtering:
    input:
        config = "{params.dirname}/bin/GFF3_extract_TranscriptID_for_filtering",
        genome_repeat = "0.RepeatMasker/genome.repeat.gff3",
        geneModels_gb = "6.combineGeneModels/geneModels.gb.gff3",
        geneModels_ge = "6.combineGeneModels/geneModels.ge.gff3",
        geneModels_gf = "6.combineGeneModels/geneModels.gf.gff3",
        filling_need_transcriptID = "filling_need_transcriptID.txt"
    output:
        protected("6.combineGeneModels/05.extract_TranscriptID_for_filtering.ok")
    params:
        dirname = "/path/to/dirname" # Replace with the correct path
    shell:
        """
        cd 6.combineGeneModels
        echo "Extracting transcript IDs for filtering..."
        {input.config} {input.genome_repeat} {input.geneModels_gb} {input.geneModels_ge} {input.geneModels_gf} > transcriptID_for_filtering.txt
        perl -ne 'print "$1\tNotEnoughEvidence\n" if m/ID=([^;]*\\.t\\d+);/' {input.geneModels_gb} >> transcriptID_for_filtering.txt
        perl -ne 'print "$1\tFilling2Uncomplete\n" if m/ID=([^;]*\\.t\\d+);/' {input.geneModels_gf} >> transcriptID_for_filtering.txt
        perl -e 'open IN, "{input.filling_need_transcriptID}"; while (<IN>) { s/.t\\d+\\n//; $gene{$_} = 1; } while (<>) { print "$1\tFilling2Complete\n" if m/ID=(([^;]+)\\.t\\d+);/ && exists $gene{$2} }' {input.geneModels_ge} >> transcriptID_for_filtering.txt
        touch 05.extract_TranscriptID_for_filtering.ok
        """
rule extract_proteins_for_filtering:
    input:
        gff3_to_protein = "{params.dirname}/bin/gff3_to_protein.pl",
        genome = "genome_path", # replace with the path to your genome file
        geneModels_gb = "6.combineGeneModels/geneModels.gb.gff3",
        geneModels_gf = "6.combineGeneModels/geneModels.gf.gff3",
        geneModels_ge = "6.combineGeneModels/geneModels.ge.gff3",
        transcriptID_for_filtering = "6.combineGeneModels/transcriptID_for_filtering.txt",
        fasta_extract_subseqs_from_list = "{params.dirname}/bin/fasta_extract_subseqs_from_list.pl"
    output:
        protected("6.combineGeneModels/06.extract_proteins_for_filtering.ok")
    params:
        dirname = "/path/to/dirname" # replace with the path to your dirname
    shell:
        """
        cd 6.combineGeneModels
        echo "Extracting proteins for filtering..."
        {input.gff3_to_protein} {input.genome} {input.geneModels_gb} {input.geneModels_gf} {input.geneModels_ge} > proteins_all.fasta 2> gff3_to_protein.log
        perl -p -i -e 's/\\*$//' proteins_all.fasta
        {input.fasta_extract_subseqs_from_list} proteins_all.fasta {input.transcriptID_for_filtering} > proteins_for_filtering.fasta 2> fasta_extract_subseqs_from_list.log
        touch 06.extract_proteins_for_filtering.ok
        """
rule validate_sequences:
    input:
        HMM_db = "path_to_HMM_db",  # replace with the path to your HMM database
        BLASTP_db = "path_to_BLASTP_db",  # replace with the path to your BLASTP database
        proteins_for_filtering = "6.combineGeneModels/proteins_for_filtering.fasta",
        para_hmmscan = "{params.dirname}/bin/para_hmmscan",
        parsing_blast_result = "{params.dirname}/bin/parsing_blast_result.pl"
    output:
        protected("7.validateSequences/07.validating.ok")
    params:
        dirname = "/path/to/dirname",  # replace with the path to your dirname
        cpu = "16"  # replace with the number of threads you want to use
    shell:
        """
        cd 7.validateSequences
        echo "Validating sequences..."
        if [ -e "{input.HMM_db}" ]; then
            for db in $(ls {input.HMM_db}*); do
                {input.para_hmmscan} --cpu {params.cpu} --hmm_db $db proteins_for_filtering.fasta >> validation_hmmscan.tab 2>> para_hmmscan.log
            done
        fi
        if [ -e "{input.BLASTP_db}" ]; then
            for db in $(ls {input.BLASTP_db}*); do
                diamond blastp --outfmt 5 --db $db --query proteins_for_filtering.fasta --out validation_blastp.xml --threads {params.cpu} &>> diamond_blastp.log
                {input.parsing_blast_result} --out-hit-confidence validation_blastp.xml >> validation_blastp.tab
            done
        else
            diamond makedb --db homolog --in ../homolog.fasta &> diamond_makedb.log
            diamond blastp --outfmt 5 --db homolog --query proteins_for_filtering.fasta --out validation_blastp.xml --threads {params.cpu} &> diamond_blastp.log
            {input.parsing_blast_result} --out-hit-confidence validation_blastp.xml > validation_blastp.tab
        fi
        touch 07.validating.ok
        """
rule get_valid_ids_and_models:
    input:
        validation_hmmscan = "7.validateSequences/validation_hmmscan.tab",
        validation_blastp = "7.validateSequences/validation_blastp.tab",
        transcriptID_for_filtering = "6.combineGeneModels/transcriptID_for_filtering.txt",
        geneModels_gb = "6.combineGeneModels/geneModels.gb.gff3",
        geneModels_ge = "6.combineGeneModels/geneModels.ge.gff3",
        geneModels_gf = "6.combineGeneModels/geneModels.gf.gff3",
        get_valid_transcriptID = "{params.dirname}/bin/get_valid_transcriptID",
        get_valid_geneModels = "{params.dirname}/bin/get_valid_geneModels"
    output:
        transcriptID_validating_passed = "8.getValidTranscriptID/transcriptID_validating_passed.tab",
        get_valid_transcriptID_ok = "8.getValidTranscriptID/08.get_valid_transcriptID.ok",
        get_valid_geneModels_ok = "9.getValidGeneModels/09.get_valid_geneModels.ok"
    params:
        dirname = "/path/to/dirname",  # replace with the path to your dirname
        get_valid_transcriptID_config = "get_valid_transcriptID_config",  # replace with your actual config file
        get_valid_geneModels_config = "get_valid_geneModels_config"  # replace with your actual config file
    shell:
        """
        cd 8.getValidTranscriptID
        echo "Getting valid transcript IDs..."
        {input.get_valid_transcriptID} {params.get_valid_transcriptID_config} {input.validation_hmmscan} {input.validation_blastp} > transcriptID_validating_passed.tab 2> get_valid_transcriptID.log
        touch 08.get_valid_transcriptID.ok

        cd ../9.getValidGeneModels
        echo "Getting valid gene models..."
        {input.get_valid_geneModels} {params.get_valid_geneModels_config} --out_prefix geneModels.h {input.transcriptID_for_filtering} {output.transcriptID_validating_passed} {input.geneModels_gb} {input.geneModels_ge} {input.geneModels_gf} 2> get_valid_geneModels.log
        touch 09.get_valid_geneModels.ok
        """
rule fillingEndsOfGeneModels:
    input:
        genome = "path/to/genome.fasta",
        geneModels = "9.getValidGeneModels/geneModels.h.coding.gff3",
        fillingEndsOfGeneModels = "{params.dirname}/bin/fillingEndsOfGeneModels"
    output:
        geneModels_filled = "10.fillingEndsOfGeneModels/geneModels.i.coding.gff3",
        fillingEndsOfGeneModels_ok = "10.fillingEndsOfGeneModels/10.fillingEndsOfGeneModels"
    params:
        dirname = "/path/to/dirname",  # replace with the path to your dirname
        fillingEndsOfGeneModels_config = "fillingEndsOfGeneModels_config"  # replace with your actual config file
    shell:
        """
        cd 10.fillingEndsOfGeneModels
        echo "Filling ends of gene models..."
        {input.fillingEndsOfGeneModels} {params.fillingEndsOfGeneModels_config} {input.genome} {input.geneModels} > geneModels.i.coding.gff3 2> fillingEndsOfGeneModels.2.log
        touch 10.fillingEndsOfGeneModels
        """

