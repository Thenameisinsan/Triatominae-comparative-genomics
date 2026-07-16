#!/usr/bin/env python3

import os
import subprocess
import sys

THREADS = 80

FASTQ_DIR = "/mnt/data2/insan/RNASeq/Fasta_files/rawReads"
FASTQC_DIR = "/mnt/data2/insan/RNASeq/Fasta_files/fastqc"
OUT_DIR = "/mnt/data2/insan/RNASeq/trimmed_reads"


def main():
    # Confirm that the input directories exist.
    for directory in [FASTQ_DIR, FASTQC_DIR]:
        if not os.path.isdir(directory):
            print(f"Directory not found: {directory}")
            sys.exit(1)

    os.makedirs(OUT_DIR, exist_ok=True)

    fastqc_folders = sorted(
        folder
        for folder in os.listdir(FASTQC_DIR)
        if folder.endswith("_1_fastqc")
        and os.path.isdir(os.path.join(FASTQC_DIR, folder))
    )

    if not fastqc_folders:
        print(f"No '_1_fastqc' directories found in {FASTQC_DIR}")
        sys.exit(1)

    print(f"Found {len(fastqc_folders)} samples.")

    for folder in fastqc_folders:
        sample = folder.removesuffix("_1_fastqc")

        print("=" * 60)
        print(f"Processing sample: {sample}")
        print("=" * 60)

        r1_fastq = os.path.join(
            FASTQ_DIR,
            f"{sample}_1.fq.gz"
        )

        r2_fastq = os.path.join(
            FASTQ_DIR,
            f"{sample}_2.fq.gz"
        )

        r1_overrep = os.path.join(
            FASTQC_DIR,
            f"{sample}_1_fastqc",
            f"{sample}_1_overrep.fasta"
        )

        r2_overrep = os.path.join(
            FASTQC_DIR,
            f"{sample}_2_fastqc",
            f"{sample}_2_overrep.fasta"
        )

        sample_out = os.path.join(OUT_DIR, sample)
        os.makedirs(sample_out, exist_ok=True)

        adapter_file = os.path.join(
            sample_out,
            f"{sample}_combined_overrep.fasta"
        )

        required_files = [
            r1_fastq,
            r2_fastq,
            r1_overrep,
            r2_overrep,
        ]

        missing_files = [
            path for path in required_files
            if not os.path.isfile(path)
        ]

        if missing_files:
            print(f"Skipping {sample}: missing input file(s).")

            for path in missing_files:
                print(f"  Missing: {path}")

            print()
            continue

        with open(adapter_file, "w") as outfile:
            for fasta_file in [r1_overrep, r2_overrep]:
                if os.path.getsize(fasta_file) == 0:
                    continue

                with open(fasta_file, "r") as infile:
                    content = infile.read().rstrip()

                    if content:
                        outfile.write(content)
                        outfile.write("\n")

        if (
            not os.path.exists(adapter_file)
            or os.path.getsize(adapter_file) == 0
        ):
            print(f"Skipping {sample}: combined adapter file is empty.\n")
            continue

        r1_paired = os.path.join(
            sample_out,
            f"{sample}_1_trimmed_paired.fq.gz"
        )

        r1_unpaired = os.path.join(
            sample_out,
            f"{sample}_1_trimmed_unpaired.fq.gz"
        )

        r2_paired = os.path.join(
            sample_out,
            f"{sample}_2_trimmed_paired.fq.gz"
        )

        r2_unpaired = os.path.join(
            sample_out,
            f"{sample}_2_trimmed_unpaired.fq.gz"
        )

        cmd = [
            "trimmomatic",
            "PE",
            "-threads",
            str(THREADS),
            "-phred33",
            r1_fastq,
            r2_fastq,
            r1_paired,
            r1_unpaired,
            r2_paired,
            r2_unpaired,
            f"ILLUMINACLIP:{adapter_file}:2:30:10",
            "SLIDINGWINDOW:4:20",
            "MINLEN:36",
        ]

        try:
            subprocess.run(cmd, check=True)
            print(f"Completed: {sample}\n")

        except FileNotFoundError:
            print(
                "Error: 'trimmomatic' was not found. "
                "Make sure it is installed and available in PATH."
            )
            sys.exit(1)

        except subprocess.CalledProcessError as error:
            print(f"Failed: {sample}")
            print(f"Return code: {error.returncode}\n")

    print("All finished.")


if __name__ == "__main__":
    main()
