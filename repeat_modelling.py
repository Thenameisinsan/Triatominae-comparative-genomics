import os
import subprocess
import glob
import sys

GENOME_DIR = "/home/insan/CompGenomics/genomes/all_genomes"
INVREP_PATH = "/home/insan/CompGenomics/genomes/all_genomes/repbase/invrep.fasta"
THREADS = "80"

BUILD_DB = "/home/insan/miniforge3/envs/repeat_tool/bin/BuildDatabase"
RMODELER = "/home/insan/miniforge3/envs/repeat_tool/bin/RepeatModeler"
RMASKER = "/home/insan/miniforge3/envs/repeat_tool/bin/RepeatMasker"


def run_pipeline():
    if not os.path.exists(INVREP_PATH):
        print(f"Can't find the invrep database at {INVREP_PATH}")
        sys.exit(1)

    search_pattern = os.path.join(GENOME_DIR, "*.fasta")
    fasta_files = glob.glob(search_pattern)

    if not fasta_files:
        print(f"No files found in {GENOME_DIR}")
        sys.exit(1)

    print(f"Found {len(fasta_files)} fasta files\n")

    for fasta in fasta_files:
        base_name = os.path.splitext(os.path.basename(fasta))[0]
        db_name = f"{base_name}_db"

        print("=" * 60)
        print(f"PROCESSING SPECIES: {base_name}")
        print("=" * 60)

        try:
            print(f"[{base_name}] Building RepeatModeler database.")
            subprocess.run(
                [BUILD_DB, "-name", db_name, fasta],
                check=True
            )

            print(f"[{base_name}] Running RepeatModeler.")
            subprocess.run(
                [
                    RMODELER,
                    "-engine", "ncbi",
                    "-threads", THREADS,
                    "-database", db_name,
                ],
                check=True
            )

            denovo_lib = f"{db_name}-families.fa"
            combined_lib = f"{base_name}_combined_library.fasta"

            print(
                f"[{base_name}] Merging de novo repeats "
                "with invrep.fasta."
            )

            if os.path.exists(denovo_lib):
                with open(combined_lib, "w") as outfile:
                    # Write the de novo repeats.
                    with open(denovo_lib, "r") as infile1:
                        outfile.write(infile1.read())
                        outfile.write("\n")

                    # Write the invertebrate repeats.
                    with open(INVREP_PATH, "r") as infile2:
                        outfile.write(infile2.read())
                        outfile.write("\n")
            else:
                print(f"{denovo_lib} not found!")
                continue

            out_dir = f"{base_name}_RM_out"

            print(f"[{base_name}] Running RepeatMasker.")
            subprocess.run(
                [
                    RMASKER,
                    "-threads", THREADS,
                    "-lib", combined_lib,
                    "-xsmall",
                    "-dir", out_dir,
                    fasta,
                ],
                check=True
            )

            print(
                f"[{base_name}] COMPLETED! "
                f"Results are in {out_dir}/\n"
            )

        except subprocess.CalledProcessError as e:
            print(f"Failed for {base_name} at command: {e.cmd}")
            print("Moving to the next genome.\n")

    print("All finished!")


if __name__ == "__main__":
    run_pipeline()
