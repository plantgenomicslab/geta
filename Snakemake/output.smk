rule output_genome:
    input:
        genome="path/to/genome.fasta"
    output:
        "path/to/output/{out_prefix}.genome.fasta"
    shell:
        """
        ln -s {input.genome} {output}
        """

rule output_gff3_gene_models:
    input:
        geneModels="6.combineGeneModels/geneModels.i.coding.gff3",
        gff3Clear="{params.dirname}/bin/GFF3Clear",
    output:
        gff3_out="path/to/output/{out_prefix}.geneModels.gff3"
    params:
        dirname="/path/to/dirname",
        out_prefix="prefix_for_output_files",
        gene_prefix="gene_prefix",
        genome="path/to/genome.fasta",
    shell:
        """
        {input.gff3Clear} --GFF3_source GETA --gene_prefix {params.gene_prefix} --gene_code_length 6 --genome {params.genome} {input.geneModels} > {output.gff3_out} 2> /dev/null
        """

rule output_best_gene_models:
    input:
        gff3_gene_models=rules.output_gff3_gene_models.output.gff3_out,
        GFF3_extract_bestGeneModels="{params.dirname}/bin/GFF3_extract_bestGeneModels",
    output:
        best_gene_models_gff3="path/to/output/{out_prefix}.bestGeneModels.gff3",
        as_num_stats="path/to/output/{out_prefix}.AS_num_of_codingTranscripts.stats",
    params:
        out_prefix="prefix_for_output_files",
    shell:
        """
        {input.GFF3_extract_bestGeneModels} {input.gff3_gene_models} > {output.best_gene_models_gff3} 2> {output.as_num_stats}
        """
rule output_gtf_gene_models:
    input:
        gff3_gene_models=rules.output_gff3_gene_models.output.gff3_out,
        gff3ToGtf="{params.dirname}/bin/gff3ToGtf.pl"
    output:
        gtf_out="path/to/output/{out_prefix}.geneModels.gtf"
    params:
        dirname="/path/to/dirname",
        out_prefix="prefix_for_output_files",
        genome="path/to/genome.fasta"
    shell:
        """
        {input.gff3ToGtf} {params.genome} {input.gff3_gene_models} > {output.gtf_out} 2> /dev/null
        """

rule output_best_gene_models_gtf:
    input:
        best_gene_models_gff3=rules.output_best_gene_models.output.best_gene_models_gff3,
        gff3ToGtf="{params.dirname}/bin/gff3ToGtf.pl"
    output:
        best_gene_models_gtf="path/to/output/{out_prefix}.bestGeneModels.gtf"
    params:
        dirname="/path/to/dirname",
        out_prefix="prefix_for_output_files",
        genome="path/to/genome.fasta"
    shell:
        """
        {input.gff3ToGtf} {params.genome} {input.best_gene_models_gff3} > {output.best_gene_models_gtf} 2> /dev/null
        """

rule output_gene_model_stats:
    input:
        best_gene_models_gtf=rules.output_best_gene_models_gtf.output.best_gene_models_gtf,
        eukaryotic_gene_model_stats="{params.dirname}/bin/eukaryotic_gene_model_statistics.pl"
    output:
        gene_models_stats="path/to/output/{out_prefix}.geneModels.stats"
    params:
        dirname="/path/to/dirname",
        out_prefix="prefix_for_output_files",
        genome="path/to/genome.fasta"
    shell:
        """
        {input.eukaryotic_gene_model_stats} {input.best_gene_models_gtf} {params.genome} {params.out_prefix} &> {output.gene_models_stats}
        """
rule output_repeat_stats:
    input:
        repeat_stats_tmp="{out_prefix}.tmp/0.RepeatMasker/genome.repeat.stats"
    output:
        repeat_stats_out="{out_prefix}.repeat.stats"
    shell:
        """
        cp {input.repeat_stats_tmp} {output.repeat_stats_out}
        """

rule output_repeat_gff3:
    input:
        repeat_gff3_tmp="{out_prefix}.tmp/0.RepeatMasker/genome.repeat.gff3"
    output:
        repeat_gff3_out="{out_prefix}.repeat.gff3"
    shell:
        """
        cp {input.repeat_gff3_tmp} {output.repeat_gff3_out}
        """

rule output_repeat_lib:
    input:
        repeat_modeler_fa="{out_prefix}.tmp/0.RepeatMasker/repeatModeler/species-families.fa",
        RM_lib="{RM_lib_path}"
    output:
        repeat_lib_out="{out_prefix}.repeat.lib"
    params:
        RM_lib_exists=os.path.isfile("{RM_lib_path}")
    run:
        if params.RM_lib_exists:
            shell("ln -sf {input.RM_lib} {output.repeat_lib_out}")
        else:
            shell("cp {input.repeat_modeler_fa} {output.repeat_lib_out}")

rule output_transcript_alignment:
    input:
        transcript_alignment_tmp="{out_prefix}.tmp/3.transcript/transfrag.alignment.gff3"
    output:
        transcript_alignment_out="{out_prefix}.transfrag_alignment.gff3"
    shell:
        """
        cp {input.transcript_alignment_tmp} {output.transcript_alignment_out}
        """

rule output_transcript_prediction:
    input:
        transcript_prediction_tmp="{out_prefix}.tmp/3.transcript/transfrag.genome.gff3"
    output:
        transcript_prediction_out="{out_prefix}.transfrag_prediction.gff3"
    shell:
        """
        cp {input.transcript_prediction_tmp} {output.transcript_prediction_out}
        """

rule output_homolog_prediction:
    input:
        homolog_prediction_tmp="{out_prefix}.tmp/4.homolog/genewise.gff3"
    output:
        homolog_prediction_out="{out_prefix}.homolog_prediction.gff3"
    shell:
        """
        cp {input.homolog_prediction_tmp} {output.homolog_prediction_out}
        """

rule output_augustus_prediction:
    input:
        augustus_prediction_tmp="{out_prefix}.tmp/5.augustus/augustus.gff3"
    output:
        augustus_prediction_out="{out_prefix}.augustus_prediction.gff3"
    params:
        dirname="{dirname}"
    shell:
        """
        {params.dirname}/bin/GFF3Clear --genome {genome} --no_attr_add {input.augustus_prediction_tmp} > {output.augustus_prediction_out}
        """
rule get_genome_repeat_stats:
    input:
        genome_repeat_stats="{out_prefix}.tmp/0.RepeatMasker/genome.repeat.stats"
    output:
        summary="{out_prefix}.gene_prediction.summary"
    shell:
        """
        awk '/^Genome Size/ || /^Repeat Ratio/ {{print $0}}' {input.genome_repeat_stats} > {output.summary}
        """

rule get_alignment_rate:
    input:
        log="{out_prefix}.tmp/2.hisat2/hisat2.log",
        previous_summary=rules.get_genome_repeat_stats.output.summary
    output:
        summary="{out_prefix}.gene_prediction.summary"
    shell:
        """
        awk '/^(\S+) overall alignment rate/ {{print "The alignment rate of RNA-Seq reads is: "$1"\n"}}' {input.log} >> {output.summary}
        """
rule get_transcript_gene_counts:
    input:
        transcript_gff3="{out_prefix}.tmp/3.transcript/transfrag.genome.gff3",
        previous_summary=rules.get_alignment_rate.output.summary
    output:
        summary="{out_prefix}.gene_prediction.summary"
    shell:
        """
        echo $(grep -c '\tgene\t' {input.transcript_gff3}) "genes were predicted by Transfrag, including" \
        $(grep -c 'Integrity=complete' {input.transcript_gff3}) "complete," \
        $(grep -c 'Integrity=5prime_partial' {input.transcript_gff3}) "5prime_partial," \
        $(grep -c 'Integrity=3prime_partial' {input.transcript_gff3}) "3prime_partial," \
        $(grep -c 'Integrity=internal' {input.transcript_gff3}) "internal genes." >> {output.summary}
        """

rule get_homolog_gene_counts:
    input:
        genewise_gff3="{out_prefix}.tmp/4.homolog/genewise.gff3",
        previous_summary=rules.get_transcript_gene_counts.output.summary
    output:
        summary="{out_prefix}.gene_prediction.summary"
    shell:
        """
        echo $(grep -c '\tgene\t' {input.genewise_gff3}) "genes were predicted by Homolog, including" \
        $(grep -c 'Integrity=complete' {input.genewise_gff3}) "complete," \
        $(grep -c 'Integrity=5prime_partial' {input.genewise_gff3}) "5prime_partial," \
        $(grep -c 'Integrity=3prime_partial' {input.genewise_gff3}) "3prime_partial," \
        $(grep -c 'Integrity=internal' {input.genewise_gff3}) "internal genes." >> {output.summary}
        """

rule check_augustus_out_exists:
    input:
        augustus_out1="{out_prefix}.tmp/5.augustus/training_again/secondtest.out",
        augustus_out2="{out_prefix}.tmp/5.augustus/training/secondtest.out",
        previous_summary=rules.get_homolog_gene_counts.output.summary
    output:
        summary="{out_prefix}.gene_prediction.summary",
        selected_augustus_out=touch("{out_prefix}.tmp/5.augustus/selected_augustus_out")
    shell:
        """
        if [ -e {input.augustus_out1} ]; then
            cp {input.augustus_out1} {output.selected_augustus_out}
        else
            cp {input.augustus_out2} {output.selected_augustus_out}
        fi
        """
rule get_augustus_accuracy_and_gene_counts:
    input:
        selected_augustus_out=rules.check_augustus_out_exists.output.selected_augustus_out,
        augustus_gff3="{out_prefix}.tmp/5.augustus/augustus.gff3",
        previous_summary=rules.check_augustus_out_exists.output.summary
    output:
        summary="{out_prefix}.gene_prediction.summary"
    shell:
        """
        ACCURACY=$(awk 'BEGIN{FS=OFS="\t"} /^nucleotide level/ {acc1=$2*3+$3*2} /^exon level/ {acc2=$2*4+$3*3} /^gene level/ {acc3=$2*2+$3*1} END{print (acc1+acc2+acc3)/15}' {input.selected_augustus_out})
        ACCURACY=$(printf "%.2f" $ACCURACY)
        echo "The accuary of AUGUSTUS Training is $ACCURACY\%." >> {output.summary}
        awk 'BEGIN{FS=OFS="\t"} /^nucleotide level/ {print "Level\tSensitivity\tSpecificity", "nucleotide_level\t"$2"\t"$3} /^exon level/ {print "exon_level\t"$2"\t"$3} /^gene level/ {print "gene_level\t"$2"\t"$3}' {input.selected_augustus_out} >> {output.summary}
        echo $(grep -c '\tgene\t' {input.augustus_gff3}) "genes were predicted by AUGUSTUS." >> {output.summary}
        """

rule get_combination_and_filter_stats:
    input:
        gene_models_a_gff3="{out_prefix}.tmp/6.combineGeneModels/geneModels.a.gff3",
        gene_models_b_gff3="{out_prefix}.tmp/6.combineGeneModels/geneModels.b.gff3",
        previous_summary=rules.get_augustus_accuracy_and_gene_counts.output.summary
    output:
        summary="{out_prefix}.gene_prediction.summary"
    shell:
        """
        echo "Statistics of the combination of 3 gene prediction methods and filtration of gene models:" >> {output.summary}
        NUM_AUGUSTUS=$(grep -c '\tgene\t' {input.gene_models_a_gff3} | grep 'augustus')
        NUM_TRANSFRAG=$(grep -c '\tgene\t' {input.gene_models_a_gff3} | grep 'transfrag')
        NUM_GENEWISE=$(grep -c '\tgene\t' {input.gene_models_a_gff3} | grep 'genewise')
        NUM_B=$(grep -c '\tgene\t' {input.gene_models_b_gff3})
        echo "(1) After first round of combination in which the AUGUSTUS results were mainly used, $NUM_AUGUSTUS genes models were supported by enough evidences, $NUM_TRANSFRAG genes models were come from transcript, $NUM_GENEWISE genes models were come from homolog, and $NUM_B genes did not supported by enough evidences." >> {output.summary}
        """

rule get_round2_combination_stats:
    input:
        picked_log="{out_prefix}.tmp/6.combineGeneModels/picked_evidence_geneModels.log",
        previous_summary=rules.get_combination_and_filter_stats.output.summary
    output:
        summary="{out_prefix}.gene_prediction.summary"
    shell:
        """
        NUMBER1=$(awk 'NR==1{print $1}' {input.picked_log})
        NUMBER2=$(awk 'NR==5{print $NF}' {input.picked_log})
        NUMBER3=$(awk 'NR==6{print $1}' {input.picked_log})
        echo "(2) In the second round of combination, $NUMBER1 evidence gene models were processed, and $NUMBER3 accurate gene models were picked out and replaced the genes predicted by AUGUSTUS, $NUMBER2 of which had the same CDS structures with the gene models predicted by AUGUSTUS." >> {output.summary}
        """

rule get_round2_combination_and_fill_stats:
    input:
        filling_log="{out_prefix}.tmp/6.combineGeneModels/fillingEndsOfGeneModels.1.log",
        previous_summary=rules.get_round2_combination_stats.output.summary
    output:
        summary="{out_prefix}.gene_prediction.summary"
    shell:
        """
        NUMBERS1=$(awk 'NR==5{print $0}' {input.filling_log})
        NUMBERS2=$(awk 'NR==6{print $0}' {input.filling_log})
        echo "(3) After the two steps of combination, ${NUMBERS1[0]} gene models with enough evidence supported were predicted. ${NUMBERS1[1]} gene models were complete; ${NUMBERS1[2]} gene models were uncomplete. ${NUMBERS2[0]} uncomplete gene models can be filled to complete, and ${NUMBERS2[1]} can not." >> {output.summary}
        """

rule get_hmm_and_blastp_validation_stats:
    input:
        valid_log="{out_prefix}.tmp/6.combineGeneModels/get_valid_geneModels.log",
        previous_summary=rules.get_round2_combination_and_fill_stats.output.summary
    output:
        summary="{out_prefix}.gene_prediction.summary"
    shell:
        """
        NUMBERS=$(awk 'NR==1{print $0}' {input.valid_log})
        echo "(4) HMM and BLASTP validation were performed to ${NUMBERS[2]} protein sequences of ${NUMBERS[1]} genes, and ${NUMBERS[5]} protein sequences of ${NUMBERS[4]} genes had valid alignment results. There are ${NUMBERS[3]} accurate gene models which did not need validation." >> {output.summary}
        """

rule get_final_gene_models:
    input:
        gene_models_gff3="{out_prefix}.geneModels.gff3",
        previous_summary=rules.get_hmm_and_blastp_validation_stats.output.summary
    output:
        summary="{out_prefix}.gene_prediction.summary"
    shell:
        """
        GENE_NUM=$(grep -c '\tgene\t' {input.gene_models_gff3})
        AS_GENE_NUM=$(grep -Po 'ID=[^;]+' {input.gene_models_gff3} | sort | uniq -c | awk '$1>=2' | wc -l)
        AUGUSTUS_NUM=$(grep -c '\tgene\t.*Source=augustus' {input.gene_models_gff3})
        TRANSFRAG_NUM=$(grep -c '\tgene\t.*Source=transfrag' {input.gene_models_gff3})
        GENEWISE_NUM=$(grep -c '\tgene\t.*Source=genewise' {input.gene_models_gff3})
        echo "(5) Finally, $GENE_NUM coding gene models were obtained, $AS_GENE_NUM of which had alternative splicing, $AUGUSTUS_NUM gene models were come from AUGUSTUS prediction, $TRANSFRAG_NUM gene models were come from transcript, $GENEWISE_NUM gene models were come from homolog." >> {output.summary}
        """

rule get_filtered_gene_models_stats:
    input:
        lncrna_gff3="{out_prefix}.geneModels_lncRNA.gff3",
        low_quality_gff3="{out_prefix}.geneModels_lowQuality.gff3",
        previous_summary=rules.get_final_gene_models.output.summary
    output:
        summary="{out_prefix}.gene_prediction.summary"
    shell:
        """
        NUM_OF_GENE1=$(grep -c '\tgene\t' {input.lncrna_gff3})
        NUM_OF_GENE2=$(grep -c '\tgene\t' {input.low_quality_gff3})
        NUM_OF_GENE3=$(($NUM_OF_GENE1 + $NUM_OF_GENE2))
        echo "(6) $NUM_OF_GENE3 gene models were filtered, and $NUM_OF_GENE1 of which had lncRNA transcripts." >> {output.summary}
        """

rule completion_message:
    input:
        previous_summary=rules.get_filtered_gene_models_stats.output.summary
    shell:
        """
        echo "\n============================================"
        echo "GETA complete successfully! ($(date))"
        echo "\n"
        """
