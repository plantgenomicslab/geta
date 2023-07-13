rule augustus_hint_preparing:
    output:
        touch("prepareAugusutusHints.ok")
    params:
        dirname="your_directory_name",  # Set this parameter accordingly
        config={
            "prepareAugusutusHints": "path/to/prepareAugusutusHints.config"
        }
    shell:
        """
        cmdString="$dirname/bin/prepareAugusutusHints ${config[prepareAugusutusHints]} ../3.transcript/intron.txt " \
                  "../3.transcript/transfrag.genome.gff3 ../4.homolog/genewise.gff3 " \
                  "../4.homolog/genewise.start_stop_hints.gff > hints.gff"
        if [ ! -e prepareAugusutusHints.ok ]; then
            echo "$(date): CMD: $cmdString"
            $cmdString
            touch prepareAugusutusHints.ok
        else
            echo "CMD(Skipped): $cmdString"
        fi
        """

rule longest_gene_length:
    output:
        "longest_gene_length.txt"
    shell:
        """
        awk '$3 == "gene" {print $4, $5}' ../3.transcript/transfrag.genome.gff3 | \
        awk '{len = $2 - $1 + 1; print len}' | sort -n | tail -n 1 > longest_gene_length.txt
        """

rule augustus_gene_prediction:
    input:
        "longest_gene_length.txt",
        "prepareAugusutusHints.ok"
    output:
        touch("first_augustus.ok")
    params:
        dirname="your_directory_name",  # Set this parameter accordingly
        config={
            "paraAugusutusWithHints": "path/to/paraAugusutusWithHints.config"
        },
        augustus_species="your_species_name",  # Set the appropriate species name
        cpu=8  # Set the desired CPU value
    shell:
        """
        gene_length=$(cat longest_gene_length.txt)
        segmentSize=5000000
        overlapSize=100000

        if ((gene_length * 4 > overlapSize)); then
            overlapSize=$((gene_length * 4))
            overlapSize_length=${#overlapSize}
            ((overlapSize_length -= 2))
            ((overlapSize_length -= 2))
            overlapSize=$(( (overlapSize / (10 ** overlapSize_length)) + 1 * (10 ** overlapSize_length) ))
            segmentSize=$((overlapSize * 50))
        fi

        cmdString="$dirname/bin/paraAugusutusWithHints ${config[paraAugusutusWithHints]} " \
                  "--species $augustus_species --cpu $cpu --segmentSize $segmentSize " \
                  "--overlapSize $overlapSize --tmp_dir aug_para_with_hints.tmp1 " \
                  "../0.RepeatMasker/genome.masked.fasta hints.gff > augustus.1.gff3"
        if [ ! -e first_augustus.ok ]; then
            echo "$(date): CMD: $cmdString"
            $cmdString
            rm -f command.augustus.list.completed
            touch first_augustus.ok
        else
            echo "CMD(Skipped): $cmdString"
        fi
        """
