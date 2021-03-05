#!/bin/bash

#PBS -W group_list=bhurwitz
#PBS -q standard
#PBS -l select=1:ncpus=1:mem=6gb
#PBS -l walltime=01:00:00
#PBS -l place=free:shared
#PBS -M aponsero@email.arizona.edu
#PBS -m bea

#positional agrument 1 = sample name
SAMPLE=$1
cwd=$(pwd)

# get results dir corresponding to the sample
VIRSORTER_DIR="$cwd/results/VirSorter/${SAMPLE}"
VIBRANT_DIR="$cwd/results/Vibrant/${SAMPLE}_renamed"

VIRSORTER="$VIRSORTER_DIR/Predicted_viral_sequences"
VIBRANT="$VIBRANT_DIR/VIBRANT_${SAMPLE}_renamed/VIBRANT_phages_${SAMPLE}_renamed"

OUT_DIR="$cwd/$2"

# get virsorter prophages names
cd $VIRSORTER

# remove potential previous files
if [[ -f "${SAMPLE}_virsorter_prophages.fasta" ]]; then
        rm "${SAMPLE}_virsorter_prophages.fasta"
fi

if [[ -f "list_prophages_VirSorter.txt" ]]; then
        rm "list_prophages_VirSorter.txt"
fi

if [[ -f "map_names_prophages.txt" ]]; then
        rm "map_names_prophages.txt"
fi

# all categories
#    cat VIRSorter_prophages_cat-4.fasta VIRSorter_prophages_cat-5.fasta VIRSorter_prophages_cat-6.fasta > ${SAMPLE}_virsorter_prophages.fasta

# cat 4 and 5 only
cat VIRSorter_prophages_cat-4.fasta VIRSorter_prophages_cat-5.fasta > ${SAMPLE}_virsorter_prophages.fasta

    grep ">" ${SAMPLE}_virsorter_prophages.fasta > list_prophages_a.txt
    sed 's/>//' list_prophages_a.txt > list_prophages.txt
    sed 's/VIRSorter_//' list_prophages.txt > list_prophages_b.txt
    sed 's/-circular.*$//' list_prophages_b.txt > list_prophages_c.txt
    sed 's/_gene_.*$//' list_prophages_c.txt > list_prophages_VirSorter.txt
    paste list_prophages.txt list_prophages_VirSorter.txt > map_names_prophages.txt
    cp list_prophages.txt $OUT_DIR/${SAMPLE}_list_virSorter.txt

    rm list_prophages.txt
    rm list_prophages_c.txt
    rm list_prophages_b.txt
    rm list_prophages_a.txt

    cp ${SAMPLE}_virsorter_prophages.fasta $OUT_DIR

# get Vibrant prophages names
    cd $VIBRANT

# remove previous run files
    if [[ -f "${SAMPLE}_vibrant_prophages.fasta" ]]; then
        rm "${SAMPLE}_vibrant_prophages.fasta"
    fi

    if [[ -f "list_prophages_Vibrant.txt" ]]; then
        rm "list_prophages_Vibrant.txt"
    fi

    if [[ -f "map_names_prophages.txt" ]]; then
        rm "map_names_prophages.txt"
    fi

    grep "_fragment_" ${SAMPLE}_renamed.phages_combined.txt > list_prophages.txt
    sed 's/_fragment_.*$//' list_prophages.txt > list_prophages_Vibrant.txt
    paste list_prophages.txt list_prophages_Vibrant.txt > map_names_prophages_o.txt
    cp list_prophages_Vibrant.txt $OUT_DIR/${SAMPLE}_list_Vibrant.txt

# merge and remove common prophages (VirSorter is selected by default)
    MERGE="$OUT_DIR/${SAMPLE}_prophage.txt"
    join -j 2 -v 1 <(sort -k 2,2 ${VIBRANT}/map_names_prophages_o.txt) <(sort -k 2,2 ${VIRSORTER}/map_names_prophages.txt) > ${VIBRANT}/map_names_prophages.txt  

    rm ${VIBRANT}/map_names_prophages_o.txt
    rm ${VIBRANT}/list_prophages.txt

    while IFS=$' ' read -r -a myArray
    do
        echo ${myArray[1]}
        grep -A 1 -w ">${myArray[1]}" ${SAMPLE}_renamed.phages_combined.fna  >> ${SAMPLE}_vibrant_prophages.fasta
    done < ${VIBRANT}/map_names_prophages.txt

    cp ${SAMPLE}_vibrant_prophages.fasta $OUT_DIR




# rename contigs in outdir
cd $OUT_DIR

f1="${SAMPLE}_virsorter_prophages.fasta"
    sed 's/VIRSorter_//' $f1 > corr1_$f1
    sed 's/_gene_.*$/_prophage/' corr1_$f1 > corr2_$f1
    sed 's/-circular//' corr2_$f1 > corr3_$f1
    rm $f1
    rm corr1_$f1
    rm corr2_$f1
    mv corr3_$f1 $f1



f2="${SAMPLE}_vibrant_prophages.fasta"
    sed 's/_fragment_.*$/_prophage/' $f2 > corr1_$f2
    rm $f2
    mv corr1_$f2 $f2


cat $f1 $f2 > ${SAMPLE}_prophages_collected.fna

awk '{
if (/>/)
 print $0 "_" ++count;
else
 print $0;
}' ${SAMPLE}_prophages_collected.fna > ${SAMPLE}_prophages_selection1.fna

rm ${SAMPLE}_prophages_collected.fna



