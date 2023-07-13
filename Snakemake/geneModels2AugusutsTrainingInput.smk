## geneModels2AugusutsTrainingInput
## https://github.com/plantgenomicslab/geta/blob/master/bin/geta.pl
## 799-881

rule augustus_gene_prediction:
    output:
        directory("5.augustus"),
        touch("5.augustus.ok")
    params:
        use_existed_augustus_species = True,  # Set this parameter accordingly
        dirname = "your_directory_name",  # Set this parameter accordingly
        cpu = 8,  # Set the desired CPU value
        genome = "path/to/genome.fasta",  # Set the path to your genome file
        config = {
            "paraCombineGeneModels": "path/to/paraCombineGeneModels.config",
            "geneModels2AugusutsTrainingInput": "path/to/geneModels2AugusutsTrainingInput.config"
        }
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

                # Prepare Augustus training input files
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

            # First Augustus HMM Training
            if [ ! -e training ]; then
                mkdir training
                [ -e training.ok ] && rm training.ok
            fi

            if [ ! -e training.ok ]; then
                cd training
                pwd=$(pwd)
                echo "PWD: $pwd"

                # Prepare Augustus training input files
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

                    # If the number of genes used for Augustus training is less than 1000, rerun geneModels2AugusutsTrainingInput
                    # with lower thresholds to increase the gene count.
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
            fi
            echo "Step 5 Done!"
        fi
        """
