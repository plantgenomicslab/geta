## geneModels2AugusutsTrainingInput
## https://github.com/plantgenomicslab/geta/blob/master/bin/geta.pl
## 799-920

rule all:
    input:
        "5.augustus.ok"

rule augustus_gene_prediction:
    output:
        directory("5.augustus"),
        touch("5.augustus.ok")
    params:
        use_existed_augustus_species=True,  # Set this parameter accordingly
        dirname="your_directory_name",  # Set this parameter accordingly
        cpu=8,  # Set the desired CPU value
        genome="path/to/genome.fasta",  # Set the path to your genome file
        config={
            "paraCombineGeneModels": "path/to/paraCombineGeneModels.config",
            "geneModels2AugusutsTrainingInput": "path/to/geneModels2AugusutsTrainingInput.config",
            "BGM2AT": "path/to/BGM2AT.config"
        },
        augustus_species_start_from="your_species_name"  # Set the appropriate species name
    shell:
        """
        echo "============================================"
        echo "Step 5: Augustus/HMM Training" "(", date, ")"
        mkdir -p 5.augustus

        if [ ! -e 5.augustus.ok ]; then
            cd 5.augustus
            pwd=$(pwd)
            echo "PWD: $pwd"

            if [ "$use_existed_augustus_species" = true ]; then
                mkdir training
                touch training.ok
                echo "Skip Augustus training for --use_existed_augustus_species parameter set"

                cd training
                pwd=$(pwd)
                echo "PWD: $pwd"

                > blank.augustus.gff3
                > blank.intron.gff

                cmdString="$dirname/bin/paraCombineGeneModels --cpu $cpu ${config[paraCombineGeneModels]} " \
                          "blank.augustus.gff3 ../../3.transcript/transfrag.genome.gff3 " \
                          "../../4.homolog/genewise.gff3 blank.intron.gff"
                echo "$(date): CMD: $cmdString"
                $cmdString

                cmdString="$dirname/bin/GFF3Clear --genome $genome combine.1.gff3 > geneModels.gff3 2> GFF3Clear.log"
                echo "$(date): CMD: $cmdString"
                $cmdString

                cd ../
            fi

            if [ ! -e training ]; then
                mkdir training
                rm -f training.ok
            fi

            if [ ! -e training.ok ]; then
                cd training
                pwd=$(pwd)
                echo "PWD: $pwd"

                > blank.augustus.gff3
                > blank.intron.gff

                cmdString="$dirname/bin/paraCombineGeneModels --cpu $cpu ${config[paraCombineGeneModels]} " \
                          "blank.augustus.gff3 ../../3.transcript/transfrag.genome.gff3 " \
                          "../../4.homolog/genewise.gff3 blank.intron.gff"
                echo "$(date): CMD: $cmdString"
                $cmdString

                cmdString="$dirname/bin/GFF3Clear --genome $genome combine.1.gff3 > geneModels.gff3 2> GFF3Clear.log"
                echo "$(date): CMD: $cmdString"
                $cmdString

                cmdString="$dirname/bin/geneModels2AugusutsTrainingInput ${config[geneModels2AugusutsTrainingInput]} " \
                          "--out_prefix ati --cpu $cpu geneModels.gff3 $genome " \
                          "1> geneModels2AugusutsTrainingInput.log 2>&1"
                if [ ! -e geneModels2AugusutsTrainingInput.ok ]; then
                    echo "$(date): CMD: $cmdString"
                    $cmdString

                    training_genes_number=$(grep -oP "Best gene Models number:\s+\K\d+" geneModels2AugusutsTrainingInput.log)
                    if [ "$training_genes_number" -lt 1000 ]; then
                        cmdString="$dirname/bin/geneModels2AugusutsTrainingInput --min_evalue 1e-9 " \
                                  "--min_identity 0.9 --min_coverage_ratio 0.9 --min_cds_num 1 " \
                                  "--min_cds_length 450 --min_cds_exon_ratio 0.40 " \
                                  "--keep_ratio_for_excluding_too_long_gene 0.99 --out_prefix ati --cpu $cpu " \
                                  "geneModels.gff3 $genome 1> geneModels2AugusutsTrainingInput.log.Loose_thresholds 2>&1"
                        echo "$(date): CMD: $cmdString"
                        $cmdString
                    fi

                    touch geneModels2AugusutsTrainingInput.ok
                else
                    echo "CMD(Skipped): $cmdString"
                fi

                > BGM2AT.log
                cmdString="$dirname/bin/BGM2AT ${config[BGM2AT]} --augustus_species_start_from $augustus_species_start_from " \
                          "--flanking_length 0 --CPU $cpu --onlytrain_GFF3 ati.filter1.gff3 ati.filter2.gff3 $genome $augustus_species " \
                          "&> BGM2AT.log"
                echo "$(date): CMD: $cmdString"
                $cmdString

                gene_info=()
                gene_length=()
                flanking_length=()
                while IFS= read -r line; do
                    if [[ $line =~ ^(\S+)\tgene\t ]]; then
                        chr="${BASH_REMATCH[1]}"
                    elif [[ $line =~ \t([-+])\t ]]; then
                        strand="${BASH_REMATCH[1]}"
                        region=()
                    elif [[ $line =~ \t(\d+)\t(\d+)\t ]]; then
                        start="${BASH_REMATCH[1]}"
                        end="${BASH_REMATCH[2]}"
                        gene_length+=("$((end - start + 1))")
                        if [ ${#region[@]} -gt 0 ]; then
                            for r in "${region[@]}"; do
                                IFS=$'\t' read -r aa bb <<<"$r"
                                gene_length+=("$((bb - aa + 1))")
                                distance=$((aa - end - 1))
                                if [ $distance -ge 50 ]; then
                                    flanking_length+=("$distance")
                                fi
                                end=$bb
                            done
                        fi
                    elif [[ $line =~ ^$ ]]; then
                        gene_info["$chr"]["$strand"]="${region[*]}"
                    else
                        region+=("$line")
                    fi
                done < geneModels.gff3

                sorted_gene_length=($(printf '%s\n' "${gene_length[@]}" | sort -n))
                sorted_flanking_length=($(printf '%s\n' "${flanking_length[@]}" | sort -n))
                mid_idx=$(( ${#sorted_flanking_length[@]} / 2 ))
                flanking_length=$(( sorted_flanking_length[mid_idx] / 8 ))
                if [ $flanking_length -ge ${sorted_gene_length[mid_idx]} ]; then
                    flanking_length=${sorted_gene_length[mid_idx]}
                fi

                cmdString="$dirname/bin/BGM2AT ${config[BGM2AT]} --augustus_species_start_from $augustus_species_start_from " \
                          "--flanking_length $flanking_length --CPU $cpu --onlytrain_GFF3 ati.filter1.gff3 ati.filter2.gff3 $genome $augustus_species " \
                          "&> BGM2AT.log"
                echo "$(date): CMD: $cmdString"
                $cmdString

                cd ../
                touch training.ok
            else
                echo "Skip Augustus training for file training.ok exists"
            fi

            echo "Step 5 Done!"
        fi
        """
