configfile: "config.yaml"

rule prepare_augustus_hints:
    version: '2.5.1'
    output:
        touch("prepareAugusutusHints.ok")
    shell:
        """
        $dirname/bin/prepareAugusutusHints {config['prepareAugusutusHints']} ../3.transcript/intron.txt ../3.transcript/transfrag.genome.gff3 ../4.homolog/genewise.gff3 ../4.homolog/genewise.start_stop_hints.gff > hints.gff
        """

rule get_longest_gene_length:
    output:
        "gene_length.txt"
    shell:
        """
        awk '/\tgene\t/ {{print $5 - $4 + 1}}' ../3.transcript/transfrag.genome.gff3 | sort -n > {output}
        """

rule determine_segment_overlap_sizes:
    input:
        "gene_length.txt"
    output:
        "segment_overlap_sizes.txt"
    shell:
        """
        gene_length=($(cat {input}))
        max_gene_length=${{gene_length[-1]}}
        segmentSize=5000000
        overlapSize=100000
        if ((max_gene_length * 4 > overlapSize)); then
            overlapSize=$((max_gene_length * 4))
            overlapSize_length=${{#overlapSize}}
            overlapSize_length=$((overlapSize_length - 2))
            overlapSize=$(((overlapSize / (10 ** overlapSize_length) + 1) * (10 ** overlapSize_length)))
            segmentSize=$((overlapSize * 50))
        fi
        echo "$segmentSize $overlapSize" > {output}
        """

rule first_augustus_gene_prediction:
    input:
        "segment_overlap_sizes.txt",
        "hints.gff"
    output:
        touch("first_augustus.ok"),
        "augustus.1.gff3"
    shell:
        """
        segment_overlap_sizes=($(cat {input.segment_overlap_sizes}))
        segmentSize=${{segment_overlap_sizes[0]}}
        overlapSize=${{segment_overlap_sizes[1]}}
        $dirname/bin/paraAugusutusWithHints {config['paraAugusutusWithHints']} --species {augustus_species} --cpu {cpu} --segmentSize $segmentSize --overlapSize $overlapSize --tmp_dir aug_para_with_hints.tmp1 ../0.RepeatMasker/genome.masked.fasta hints.gff > augustus.1.gff3
        """

rule filter_genes_for_training:
    input:
        "augustus.1.gff3"
    output:
        "geneModels.gff3"
    shell:
        """
        awk -v keep_num=0 -v total_num=0 '
            BEGIN {
                FS = "\t"
                OFS = "\t"
                keep = 0
            }
            $3 == "gene" {
                total_num++
                gene_id = gensub(/.*ID=([^;]+);.*/, "\\1", "g", $9)
                gene_models[gene_id] = $0
                if ($9 ~ /exonHintRatio=[0-9.]+/ && substr($9, RSTART, RLENGTH) >= 95) {
                    keep = 1
                    keep_num++
                } else if ($9 ~ /intronSupport=(\\d+)\\/\\1/ && substr($9, RSTART, RLENGTH) == $9 && $9 > 0) {
                    keep = 1
                    keep_num++
                } else {
                    keep = 0
                }
                if (keep == 1) {
                    keep_genes[++num_keep_genes] = gene_id
                }
            }
            $3 == "mRNA" {
                gene_id = gensub(/.*Parent=([^;]+);.*/, "\\1", "g", $9)
                mRNA_id = gensub(/.*ID=([^;]+);.*/, "\\1", "g", $9)
                gene_models_mRNA[gene_id][mRNA_id] = $0
                mRNA_score[gene_id][mRNA_id] = $6
            }
            END {
                print "Total genes predicted by Augustus: " total_num
                print "Good genes picked for next training: " keep_num
                for (i in keep_genes) {
                    gene_id = keep_genes[i]
                    print gene_models[gene_id]
                    max_mRNA_id = ""
                    max_score = 0
                    for (mRNA_id in gene_models_mRNA[gene_id]) {
                        if (mRNA_score[gene_id][mRNA_id] > max_score) {
                            max_score = mRNA_score[gene_id][mRNA_id]
                            max_mRNA_id = mRNA_id
                        }
                    }
                    print gene_models_mRNA[gene_id][max_mRNA_id]
                    print ""
                }
            }
        ' {input} > {output}
        """

rule train_augustus_iteration:
    input:
        "geneModels.gff3",
        "../training/secondtest.out"
    output:
        touch("training_again.ok")
    params:
        enable_augustus_training_iteration=True,
        use_existed_augustus_species=False,
        dirname="your_directory_name",  # Set this parameter accordingly
        config={
            "geneModels2AugusutsTrainingInput": "path/to/geneModels2AugusutsTrainingInput.config",
            "BGM2AT": "path/to/BGM2AT.config",
            "paraAugusutusWithHints": "path/to/paraAugusutusWithHints.config"
        },
        cpu=8  # Set the desired CPU value
    shell:
        """
        if [ ! -e training_again ]; then
            mkdir training_again
            species_config_dir=$(echo $AUGUSTUS_CONFIG_PATH)
            species_config_dir="${species_config_dir}/species/${augustus_species}"
            cmdString="rm -rf $species_config_dir && cp -a ./training/hmm_files_bak/ $species_config_dir"
            echo "CMD: $cmdString"
            $cmdString
            rm -f training_again.ok
        fi
        if [ ! -e training_again.ok ]; then
            cd training_again
            pwd=$(pwd)
            echo "PWD: $pwd"

            awk -v keep_num=0 -v total_num=0 '
                BEGIN {
                    FS = "\t"
                    OFS = "\t"
                    keep = 0
                }
                $3 == "gene" {
                    total_num++
                    gene_id = gensub(/.*ID=([^;]+);.*/, "\\1", "g", $9)
                    gene_models[gene_id] = $0
                    if ($9 ~ /exonHintRatio=[0-9.]+/ && substr($9, RSTART, RLENGTH) >= 95) {
                        keep = 1
                        keep_num++
                    } else if ($9 ~ /intronSupport=(\\d+)\\/\\1/ && substr($9, RSTART, RLENGTH) == $9 && $9 > 0) {
                        keep = 1
                        keep_num++
                    } else {
                        keep = 0
                    }
                    if (keep == 1) {
                        keep_genes[++num_keep_genes] = gene_id
                    }
                }
                $3 == "mRNA" {
                    gene_id = gensub(/.*Parent=([^;]+);.*/, "\\1", "g", $9)
                    mRNA_id = gensub(/.*ID=([^;]+);.*/, "\\1", "g", $9)
                    gene_models_mRNA[gene_id][mRNA_id] = $0
                    mRNA_score[gene_id][mRNA_id] = $6
                }
                END {
                    print "Total genes predicted by Augustus: " total_num
                    print "Good genes picked for next training: " keep_num
                    for (i in keep_genes) {
                        gene_id = keep_genes[i]
                        print gene_models[gene_id]
                        max_mRNA_id = ""
                        max_score = 0
                        for (mRNA_id in gene_models_mRNA[gene_id]) {
                            if (mRNA_score[gene_id][mRNA_id] > max_score) {
                                max_score = mRNA_score[gene_id][mRNA_id]
                                max_mRNA_id = mRNA_id
                            }
                        }
                        print gene_models_mRNA[gene_id][max_mRNA_id]
                        print ""
                    }
                }
            ' ../augustus.1.gff3 > geneModels.gff3
            echo "Total genes predicted by Augustus: $(awk 'BEGIN{FS="\t";OFS="\t";total=0} $3=="gene"{total++}END{print total}' geneModels.gff3)"
            echo "Good genes picked for next traning: $(grep -c "^gene" geneModels.gff3)"

            cmdString="$dirname/bin/geneModels2AugusutsTrainingInput --out_prefix ati --cpu $cpu ${config[geneModels2AugusutsTrainingInput]} geneModels.gff3 $genome > geneModels2AugusutsTrainingInput.log 2>&1"
            if [ ! -e geneModels2AugusutsTrainingInput.ok ]; then
                echo "$(date): CMD: $cmdString"
                $cmdString
                training_genes_number=$(awk '/Best gene Models number:/{print $NF}' geneModels2AugusutsTrainingInput.log)
                if ((training_genes_number < 1000)); then
                    cmdString="$dirname/bin/geneModels2AugusutsTrainingInput --min_evalue 1e-9 --min_identity 0.9 --min_coverage_ratio 0.9 --min_cds_num 1 --min_cds_length 450 --min_cds_exon_ratio 0.40 --keep_ratio_for_excluding_too_long_gene 0.99 --out_prefix ati --cpu $cpu geneModels.gff3 $genome > geneModels2AugusutsTrainingInput.log.Loose_thresholds 2>&1"
                    echo "$(date): CMD: $cmdString"
                    $cmdString
                fi
                touch geneModels2AugusutsTrainingInput.ok
            else
                echo "CMD(Skipped): $cmdString"
            fi

            # Augustus training
            declare -A gene_info
            declare -a flanking_length
            declare -a gene_length

            while read -r line; do
                if [[ $line =~ ^gene ]]; then
                    IFS=$'\t' read -r -a fields <<< "$line"
                    chr="${fields[0]}"
                    strand="${fields[6]}"
                    region="${fields[3]}\t${fields[4]}"
                    gene_info[$chr][$strand][$region]=1
                fi
            done < geneModels.gff3

            for chr in "${!gene_info[@]}"; do
                for strand in "${!gene_info[$chr][@]}"; do
                    regions=()
                    while IFS= read -r region; do
                        regions+=("$region")
                    done < <(printf '%s\n' "${!gene_info[$chr][$strand][@]}" | sort -n)
                    first_region="${regions[0]}"
                    IFS=$'\t' read -r -a start_end <<< "$first_region"
                    start="${start_end[0]}"
                    end="${start_end[1]}"
                    gene_length+=("$((end - start + 1))")
                    for ((i = 1; i < ${#regions[@]}; i++)); do
                        IFS=$'\t' read -r -a aa_bb <<< "${regions[i]}"
                        aa="${aa_bb[0]}"
                        bb="${aa_bb[1]}"
                        gene_length+=("$((bb - aa + 1))")
                        distance=$((aa - end - 1))
                        if ((distance >= 50)); then
                            flanking_length+=("$distance")
                        fi
                        start="$aa"
                        end="$bb"
                    done
                done
            done

            IFS=$'\n' read -rd '' -a gene_length <<< "$(printf '%s\n' "${gene_length[@]}" | sort -n)"
            IFS=$'\n' read -rd '' -a flanking_length <<< "$(printf '%s\n' "${flanking_length[@]}" | sort -n)"
            flanking_length=$((flanking_length[${#flanking_length[@]} / 2] / 8))
            if ((flanking_length >= gene_length[${#gene_length[@]} / 2])); then
                flanking_length="${gene_length[${#gene_length[@]} / 2]}"
            fi

            cmdString="$dirname/bin/BGM2AT ${config[BGM2AT]} --flanking_length $flanking_length --CPU $cpu --onlytrain_GFF3 ati.filter1.gff3 --stopAfterFirstEtraining ati.filter2.gff3 $genome $augustus_species > BGM2AT.log 2>&1"
            if [ ! -e firsttest.out ]; then
                echo "$(date): CMD: $cmdString"
                $cmdString
            else
                echo "CMD(Skipped): $cmdString"
            fi

            first_accuracy=0
            awk '/nucleotide level/ || /exon level/ || /gene level/ {
                sum = ($NF * 3) + ($(NF-1) * 2)
                if (/exon level/ || /gene level/) {
                    sum += ($NF * 4) + ($(NF-1) * 3)
                }
                if (/gene level/) {
                    sum += ($NF * 2) + ($(NF-1) * 1)
                }
                first_accuracy += sum
            }
            END {
                print "The accuracy value of augustus training is: " first_accuracy / 15
            }' ../training/secondtest.out

            second_accuracy=0
            awk '/nucleotide level/ || /exon level/ || /gene level/ {
                sum = ($NF * 3) + ($(NF-1) * 2)
                if (/exon level/ || /gene level/) {
                    sum += ($NF * 4) + ($(NF-1) * 3)
                }
                if (/gene level/) {
                    sum += ($NF * 2) + ($(NF-1) * 1)
                }
                second_accuracy += sum
            }
            END {
                print "The accuracy value of augustus training iteration is: " second_accuracy / 15
            }' firsttest.out

            if ((second_accuracy > first_accuracy)); then
                cmdString="$dirname/bin/BGM2AT ${config[BGM2AT]} --flanking_length $flanking_length --CPU $cpu --onlytrain_GFF3 ati.filter1.gff3 ati.filter2.gff3 $genome $augustus_species > BGM2AT.log 2>&1"
                echo "$(date): CMD: $cmdString"
                $cmdString

                cd ..
                pwd=$(pwd)
                echo "PWD: $pwd"
                touch training_again.ok

                # Get the longest gene length
                longest_gene_length=$(awk '/\tgene\t/ {len = $5 - $4 + 1} END {print len}' augustus.1.gff3)

                segmentSize=5000000
                overlapSize=100000
                if ((longest_gene_length * 4 > overlapSize)); then
                    overlapSize=$((longest_gene_length * 4))
                    overlapSize_length=${#overlapSize}
                    ((overlapSize_length -= 2))
                    ((overlapSize_length -= 2))
                    overlapSize=$(( (overlapSize / (10 ** overlapSize_length)) + 1 * (10 ** overlapSize_length) ))
                    segmentSize=$((overlapSize * 50))
                fi

                cmdString="$dirname/bin/paraAugusutusWithHints ${config[paraAugusutusWithHints]} --species $augustus_species --cpu $cpu --segmentSize $segmentSize --overlapSize $overlapSize --tmp_dir aug_para_with_hints.tmp2 ../0.RepeatMasker/genome.masked.fasta hints.gff > augustus.2.gff3"
                echo "$(date): CMD: $cmdString"
                $cmdString

                cmdString="$dirname/bin/addHintRatioToAugustusResult training/geneModels.gff3 hints.gff augustus.2.gff3 > augustus.gff3"
                echo "$(date): CMD: $cmdString"
                $cmdString
            else
                echo "The iteration step could not increase the accuracy of augustus training! and the hmm files will be rolled back!"
                cd ..
                pwd=$(pwd)
                echo "PWD: $pwd"
                species_config_dir=$(echo $AUGUSTUS_CONFIG_PATH)
                species_config_dir="${species_config_dir}/species/${augustus_species}"
                cmdString="rm -rf $species_config_dir && cp -a ./training/hmm_files_bak/ $species_config_dir"
                echo "CMD: $cmdString"
                $cmdString
                touch training_again.ok

                cmdString="$dirname/bin/addHintRatioToAugustusResult training/geneModels.gff3 hints.gff augustus.1.gff3 > augustus.gff3"
                echo "$(date): CMD: $cmdString"
                $cmdString
            fi
        else
            echo "Skip Augustus training again for file training_again.ok exists"
            if [ ! -e augustus.gff3 ]; then
                # Get the longest gene length
                longest_gene_length=$(awk '/\tgene\t/ {len = $5 - $4 + 1} END {print len}' augustus.1.gff3)

                segmentSize=5000000
                overlapSize=100000
                if ((longest_gene_length * 4 > overlapSize)); then
                    overlapSize=$((longest_gene_length * 4))
                    overlapSize_length=${#overlapSize}
                    ((overlapSize_length -= 2))
                    ((overlapSize_length -= 2))
                    overlapSize=$(( (overlapSize / (10 ** overlapSize_length)) + 1 * (10 ** overlapSize_length) ))
                    segmentSize=$((overlapSize * 50))
                fi

                cmdString="$dirname/bin/paraAugusutusWithHints ${config[paraAugusutusWithHints]} --species $augustus_species --cpu $cpu --segmentSize $segmentSize --overlapSize $overlapSize --tmp_dir aug_para_with_hints.tmp2 ../0.RepeatMasker/genome.masked.fasta hints.gff > augustus.2.gff3"
                echo "$(date): CMD: $cmdString"
                $cmdString

                cmdString="$dirname/bin/addHintRatioToAugustusResult training/geneModels.gff3 hints.gff augustus.2.gff3 > augustus.gff3"
                echo "$(date): CMD: $cmdString"
                $cmdString
            fi
        fi
        cd ..
        touch 5.augustus.ok
        """

rule final_augustus_gene_prediction:
    input:
        "5.augustus.ok"
    output:
        "augustus.gff3"
    shell:
        """
        cp training_again/augustus.gff3 {output}
        """

# Define the workflow
rule augustus_training_iteration:
    input:
        "prepareAugusutusHints.ok",
        "first_augustus.ok",
        "training_again.ok",
        "5.augustus.ok"
    output:
        touch("augustus_training_iteration.ok")
    shell:
        """
        # Placeholder rule to signify the completion of the Augustus training iteration
        """
