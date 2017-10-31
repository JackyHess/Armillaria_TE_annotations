#!/bin/bash


#make windowed file
for filename in *.genome; do
	bedtools makewindows -w 10000 -g $filename > $(basename "$filename" .genome)_10k.bed
	bedtools makewindows -w 50000 -g $filename > $(basename "$filename" .genome)_50k.bed
	bedtools makewindows -w 100000 -g $filename > $(basename "$filename" .genome)_100k.bed
done

# remove genes overlapping with TEs more than 20%
for filename in *.genome; do
	species=$(basename "$filename" .genome)
	bedtools intersect -a ${species}_genes.bed -b ${species}_TEanno.bed -f 0.2 -v > ${species}_genes_noTE.bed
	echo $species
	wc -l ${species}_genes.bed
	wc -l ${species}_genes_noTE.bed
done



# get coverage per window
for filename in *.genome; do
	species=$(basename "$filename" .genome)
	
	# TE coverage
	bedtools coverage -a ${species}_10k.bed -b ${species}_TEanno.bed > ${species}_TEcov_10k.bed
	bedtools coverage -a ${species}_50k.bed -b ${species}_TEanno.bed > ${species}_TEcov_50k.bed
	bedtools coverage -a ${species}_100k.bed -b ${species}_TEanno.bed > ${species}_TEcov_100k.bed
	
	# Gene coverage
	bedtools coverage -a ${species}_10k.bed -b ${species}_genes_noTE.bed > ${species}_genecov_10k.bed
	bedtools coverage -a ${species}_50k.bed -b ${species}_genes_noTE.bed > ${species}_genecov_50k.bed
	bedtools coverage -a ${species}_100k.bed -b ${species}_genes_noTE.bed > ${species}_genecov_100k.bed
	
done

### combining and plotting data with R (see corresponding scripts)
### continue analysis with 50k windows since they seem to provide an informative view

# get % TE bases in TE-only windows, % genes in gene-only windows
for filename in *.genome; do
    species=$(basename "$filename" .genome)

    #TE-only coverage
    awk 'NR > 1 && $4 > 0.0 && $5 == 0.0' ${species}_mergedCov_50k.tab > ${species}_TEonly_50k.bed
    bedtools intersect -a ${species}_TEanno.bed -b ${species}_TEonly_50k.bed > ${species}_TEanno_TEonly50k.bed
    echo $species
    awk '{a = ($3 - $2)+1; b+=a} END {print "All_TE " b}' ${species}_TEanno.bed
    awk '{a = ($3 - $2)+1; b+=a} END {print "TEonly " b}' ${species}_TEanno_TEonly50k.bed

    # Gene numbers in genic regions
    awk 'NR > 1 && $5 > 0.0 && $4 == 0.0' ${species}_mergedCov_50k.tab > ${species}_genesonly_50k.bed
    bedtools intersect -a ${species}_genes_noTE.bed -b ${species}_genesonly_50k.bed > ${species}_genes_genesonly_50k.bed
    cat ${species}_genes_noTE.bed | cut -f 4 | uniq | wc -l
    cat ${species}_genes_genesonly_50k.bed | cut -f 4 | uniq | wc -l

done

# produce gene lists for functional enrichment test. Cutoffs (TE proportions): 0.2, 0.5, 0.8
for filename in *.genome; do
    species=$(basename "$filename" .genome)

    echo $species

    # At least 20% TE regions
    awk 'NR > 1 && $4 >= 0.2' ${species}_mergedCov_50k.tab > ${species}_TEgt0.2_50k.bed
    bedtools intersect -a ${species}_genes_noTE.bed -b ${species}_TEgt0.2_50k.bed | cut -f 4 | uniq > ${species}_TEgt0.2_gene.ids
    wc -l ${species}_TEgt0.2_gene.ids


    # At least 50% TE regions
    awk 'NR > 1 && $4 >= 0.5' ${species}_mergedCov_50k.tab > ${species}_TEgt0.5_50k.bed
    bedtools intersect -a ${species}_genes_noTE.bed -b ${species}_TEgt0.5_50k.bed | cut -f 4 | uniq > ${species}_TEgt0.5_gene.ids
    wc -l ${species}_TEgt0.5_gene.ids

    # At least 80% TE regions
    awk 'NR > 1 && $4 >= 0.8' ${species}_mergedCov_50k.tab > ${species}_TEgt0.8_50k.bed
    bedtools intersect -a ${species}_genes_noTE.bed -b ${species}_TEgt0.8_50k.bed | cut -f 4 | uniq > ${species}_TEgt0.8_gene.ids
    wc -l ${species}_TEgt0.8_gene.ids
done


