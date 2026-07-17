#!/bin/bash

FASTA="nelson_anoph_p450.fasta"
OUTDIR="clan_hmms"
mkdir -p $OUTDIR

declare -A CLANS

CLANS["CYP6_detox"]="^>CYPf|^>CYPh|^>CYPm|^>CYPq|^>CYPr|^>CYPp"
CLANS["CYP4_lipid"]="^>CYPa|^>CYPb|^>CYPc|^>CYPd|^>CYPe|^>CYPi|^>CYPk|^>CYPn|^>CYPo"
CLANS["CYP2_hormone"]="^>CYPs"
CLANS["CYP12_mitochondrial"]="^>CYPl"
CLANS["CYP305_JH"]="^>CYPj"
CLANS["CYP304_pheromone"]="^>CYPg"

for clan in "${!CLANS[@]}"; do
    pattern="${CLANS[$clan]}"

    echo "Processing clan: ${clan}"

    # Extract matching sequence IDs
    grep -E "$pattern" $FASTA | \
        sed 's/>//' | \
        awk '{print $1}' > ${OUTDIR}/${clan}_ids.txt

    count=$(wc -l < ${OUTDIR}/${clan}_ids.txt)
    echo "  Sequences found: ${count}"

    if [ "$count" -lt 2 ]; then
        continue
    fi

    seqtk subseq $FASTA ${OUTDIR}/${clan}_ids.txt \
        > ${OUTDIR}/${clan}.fasta

    echo "  Running MAFFT alignment..."
    mafft --auto \
        --quiet \
        --thread 8 \
        ${OUTDIR}/${clan}.fasta \
        > ${OUTDIR}/${clan}_aln.fasta

    echo "  Building HMM..."
    hmmbuild --cpu 8 \
        -n ${clan} \
        ${OUTDIR}/${clan}.hmm \
        ${OUTDIR}/${clan}_aln.fasta

done

echo "Combining clan HMMs into single database"

cat ${OUTDIR}/CYP6_detox.hmm \
    ${OUTDIR}/CYP4_lipid.hmm \
    ${OUTDIR}/CYP2_hormone.hmm \
    ${OUTDIR}/CYP12_mitochondrial.hmm \
    ${OUTDIR}/CYP305_JH.hmm \
    ${OUTDIR}/CYP304_pheromone.hmm \
    > nelson_cyp_clans_combined.hmm

hmmpress nelson_cyp_clans_combined.hmm

echo "Scanning P450 OG sequences"

mkdir -p cyp_clan_results

for OG in OG0000064 OG0000159 OG0000262 OG0000525 OG0001130 OG0001725; do
    echo "Scanning ${OG}"

    hmmscan --tblout cyp_clan_results/${OG}_clans.txt \
        --noali \
        --cpu 8 \
        -E 1e-3 \
        nelson_cyp_clans_combined.hmm \
        /home/insan/CompGenomics/miscellaneous/p450_test/fasta/${OG}.fa
done


for OG in OG0000064 OG0000159 OG0000262 OG0000525 OG0001130 OG0001725; do
    echo "${OG}"
    grep -v "^#" cyp_clan_results/${OG}_clans.txt | \
        awk '{print $1, $2, $5, $6}' | \
        sort -k3 -g | \
        awk '!seen[$1]++ {printf "%-22s %-22s %-10s %s\n", $2, $1, $3, $4}'
done


for OG in OG0000064 OG0000159 OG0000262 OG0000525 OG0001130 OG0001725; do
    echo "${OG}"
    grep -v "^#" cyp_clan_results/${OG}_clans.txt | \
        awk '!seen[$1]++ {print $2}' | \
        sort | uniq -c | sort -rn
done

echo "Done."
