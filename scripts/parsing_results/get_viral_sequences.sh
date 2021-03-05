#!/bin/bash

#PBS -W group_list=bhurwitz
#PBS -q standard
#PBS -l select=1:ncpus=2:mem=10gb
#PBS -l walltime=24:00:00
#PBS -M aponsero@email.arizona.edu
#PBS -m bea
#PBS -l place=pack:shared

cwd=$(pwd)

SAMPLE=$1
OUT_DIR="$cwd/$2"


##############################################
# virSorter parsing

    VIRSORTER_DIR="$cwd/results/VirSorter/${SAMPLE}"
    VIRSORTER="$VIRSORTER_DIR/Predicted_viral_sequences"

    cd $VIRSORTER

# remove potential previous files
    if [[ -f "${SAMPLE}_virsorter.fasta" ]]; then
        rm "${SAMPLE}_virsorter.fasta"
    fi

    if [[ -f "${SAMPLE}_list_VirSorter.txt" ]]; then
        rm "list_VirSorter.txt"
    fi

    if [[ -f "map_names.txt" ]]; then
        rm "map_names.txt"
    fi

# get viral sequences names
    cat VIRSorter_cat-1.fasta VIRSorter_cat-2.fasta VIRSorter_cat-3.fasta > ${SAMPLE}_virsorter.fasta

    grep ">" ${SAMPLE}_virsorter.fasta > list_a.txt
    sed 's/>//' list_a.txt > list.txt
    sed 's/VIRSorter_//' list.txt > list_b.txt
    sed 's/-circular.*$//' list_b.txt > list_c.txt
    sed 's/-cat_.*$//' list_c.txt > ${SAMPLE}_list_VirSorter.txt
    paste list.txt ${SAMPLE}_list_VirSorter.txt > map_names.txt

    rm list.txt
    rm list_c.txt
    rm list_b.txt
    rm list_a.txt

    cp ${SAMPLE}_list_VirSorter.txt $OUT_DIR


###############################################################
# vibrant parsing

# get Vibrant prophages names
    VIBRANT_DIR="$cwd/results/Vibrant/${SAMPLE}_renamed"
    VIBRANT="$VIBRANT_DIR/VIBRANT_${SAMPLE}_renamed/VIBRANT_phages_${SAMPLE}_renamed"

    cd $VIBRANT

# remove previous run files
    if [[ -f "${SAMPLE}_vibrant.fasta" ]]; then
        rm "${SAMPLE}_vibrant.fasta"
    fi

    if [[ -f "${SAMPLE}_list_Vibrant.txt" ]]; then
        rm "${SAMPLE}_list_Vibrant.txt"
    fi

    grep -v "_fragment_" ${SAMPLE}_renamed.phages_combined.txt > ${SAMPLE}_list_Vibrant.txt 

    cp ${SAMPLE}_list_Vibrant.txt $OUT_DIR


###############################################################
# MARVEL parsing

# get MARVEL names
    MARVEL_DIR="$cwd/results/marvel"
    MARVEL="$MARVEL_DIR/${SAMPLE}_bin/results/phage_bins"

    if [[ ! -d "$MARVEL_DIR/${SAMPLE}_bin/results/phage_bins" ]]; then
        echo "no phages found for ${SAMPLE}"

    else

        cd $MARVEL

# remove previous run files
        if [[ -f "${SAMPLE}_marvel.fasta" ]]; then
            rm "${SAMPLE}_marvel.fasta"
        fi

        if [[ -f "${SAMPLE}_list_marvel.txt" ]]; then
            rm "${SAMPLE}_list_marvel.txt"
        fi

        cat bin.* > ${SAMPLE}_marvel.fasta 
        grep ">" ${SAMPLE}_marvel.fasta > list_marvel.txt
        sed 's/>//' list_marvel.txt > ${SAMPLE}_list_marvel.txt

        rm list_marvel.txt 

        cp ${SAMPLE}_list_marvel.txt $OUT_DIR
    fi

###############################################################
# VirFinder parsing

# get VirFinder names
    VIRFINDER_DIR="$cwd/results/VirFinder"

    cd $VIRFINDER_DIR

# remove previous run files

    if [[ -f "${SAMPLE}_list_virFinder.txt" ]]; then
        rm "${SAMPLE}_list_virFinder.txt"
    fi

    awk '{ if (($4 >= 0.960) && ($3 >= 1000)) { print $2 } }' ${SAMPLE}.txt | sed 's/"//g' > ${SAMPLE}_list_virFinder.txt

    cp ${SAMPLE}_list_virFinder.txt $OUT_DIR


###############################################################
# merge viral sequences

cd $OUT_DIR
    cat ${SAMPLE}_list_* > ${SAMPLE}_merge.txt
    sort ${SAMPLE}_merge.txt | uniq > ${SAMPLE}_unique.txt 


###############################################################
# retrieve lytic sequences

CONT_DIR="$cwd/test"
cd $CONT_DIR
out="$OUT_DIR/${SAMPLE}_phages_selection1.fna"
file_id="$CONT_DIR/${SAMPLE}_renamed.fasta"

# convert multiline into singleline
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $file_id > ${file_id}_singleline.fasta

while read -r contig_id
do 
    echo "$contig_id in $file_id"
    file="${file_id}_singleline.fasta"
    grep -A1 -w $contig_id $file >> $out
done < $OUT_DIR/${SAMPLE}_unique.txt








