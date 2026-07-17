import os
import subprocess

fa_dir = '/home/insan/CompGenomics/cafe/cafe_significant_annotations/OGs'
out_path = '/home/insan/CompGenomics/cafe/cafe_significant_annotations/interpro'
ipr_script = '/home/insan/my_interproscan/interproscan-5.76-107.0/interproscan.sh'
error_log = os.path.join(out_path, "error.log")
extension = ".fa"


os.makedirs(out_path, exist_ok=True)

print(f"Scanning directory for FASTA files: {fa_dir}")

all_files = [f for f in os.listdir(fa_dir) if f.endswith(extension)]

if not all_files:
    print(f"No files found in {fa_dir} ending with {extension}")
else:
    print(f"Found {len(all_files)} files to process.")

for filename in all_files:
    og_id = os.path.splitext(filename)[0]

    input_fasta = os.path.join(fa_dir, filename)
    current_out_dir = os.path.join(out_path, og_id)
    done_file = os.path.join(current_out_dir, "done")

    os.makedirs(current_out_dir, exist_ok=True)

    if os.path.exists(done_file):
        print(f"Skipping {og_id}: Already done.")
        continue

    print(f"Running InterProScan for {og_id}...")

    output_base = os.path.join(current_out_dir, og_id)

    cmd = [
        ipr_script,
        "-i", input_fasta,
        "-appl", "Pfam,PANTHER,CDD,SMART,SUPERFAMILY",
        "-f", "tsv",
        "-goterms",
        "-pa",
        "-iprlookup",
        "-cpu", "80",
        "-b", output_base
    ]

    try:
        subprocess.run(cmd, check=True)

        tsv_result = os.path.join(current_out_dir, f"{og_id}.tsv")

        if os.path.exists(tsv_result):
            with open(done_file, "w") as f:
                f.write("done\n")

            print(f"Success: {og_id}")

        else:
            print(f"InterProScan finished but TSV not found for {og_id}")
            with open(error_log, 'a') as logfile:
                logfile.write(f"TSV missing after run: {og_id}\n")

    except subprocess.CalledProcessError as e:
        print(f"Failed: {og_id}")
        with open(error_log, 'a') as logfile:
            logfile.write(f"Failed run {og_id}: {e}\n")

print("\nAll processing finished.")
