#!/usr/bin/env python3

import os
import subprocess

THREADS = 80
TRIM_DIR = "/mnt/data2/insan/RNASeq/trimmed_reads/trimmed_reads"
OUT_DIR = "/mnt/data2/insan/RNASeq/assembly/"

species_ls = [
    ("rhodnius", "/home/insan/CompGenomics/genomes/genomes/rhodnius.genome.fasta", "rhodnius"),
    ("rubida", "/home/insan/CompGenomics/genomes/genomes/rubida.genome.fasta", "rubida")
]

os.makedirs(OUT_DIR, exist_ok=True)

for species, genome, prefix in species_ls:

    print(f"\nProcessing {species}")

    index = f"{OUT_DIR}/{species}_index"
    bam_dir = f"{OUT_DIR}/{species}_bam"
    trinity_out = f"{OUT_DIR}/{species}_trinity"

    os.makedirs(bam_dir, exist_ok=True)

    subprocess.run(f"hisat2-build {genome} {index}", shell=True, check=True)

    bams = []

    for sample in sorted(os.listdir(TRIM_DIR)):

        if not sample.startswith(prefix):
            continue

        r1 = f"{TRIM_DIR}/{sample}/{sample}_1_trimmed_paired.fq.gz"
        r2 = f"{TRIM_DIR}/{sample}/{sample}_2_trimmed_paired.fq.gz"
        bam = f"{bam_dir}/{sample}.sorted.bam"

        print(f"Mapping {sample}")

        subprocess.run(
            f"hisat2 -p {THREADS} --dta -x {index} -1 {r1} -2 {r2} | "
            f"samtools sort -@ {THREADS} -o {bam}",
            shell=True,
            check=True
        )

        subprocess.run(f"samtools index {bam}", shell=True, check=True)
        bams.append(bam)

    merged_bam = f"{bam_dir}/{species}.merged.sorted.bam"

    subprocess.run(
        f"samtools merge -f -@ {THREADS} {merged_bam} " + " ".join(bams),
        shell=True,
        check=True
    )

    subprocess.run(f"samtools index {merged_bam}", shell=True, check=True)

    subprocess.run(
        f"Trinity --genome_guided_bam {merged_bam} "
        f"--genome_guided_max_intron 10000 "
        f"--max_memory 100G "
        f"--CPU {THREADS} "
        f"--output {trinity_out}",
        shell=True,
        check=True
    )

print("\nAll finished")
