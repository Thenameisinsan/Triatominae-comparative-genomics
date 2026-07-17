import os
import subprocess

fa_dir = '/home/insan/CompGenomics/cafe/cafe_significant_annotations/OGs'
out_path = '/home/insan/CompGenomics/cafe/cafe_significant_annotations/emapper'
emap_env_path = '/home/insan/miniconda3/envs/eggnog'
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

    current_out_dir = os.path.join(out_path, og_id)
    input_fasta = os.path.join(fa_dir, filename)
    done_file = os.path.join(current_out_dir, "done")

    os.makedirs(current_out_dir, exist_ok=True)

    if os.path.exists(done_file):
        print(f"Skipping {og_id}: Already done.")
        continue

    print(f"Running emapper for {og_id}...")

    cmd = (
        f"conda run -p {emap_env_path} emapper.py"
        f"-i {input_fasta}"
        f"--itype proteins"
        f"-m diamond"
        f"--cpu 80"
        f"--output {og_id}"
        f"--output_dir {current_out_dir}"
        f"&& touch {done_file}"
    )

    try:
        subprocess.run(cmd, shell=True, check=True)
        print(f"Success: {og_id}")

    except subprocess.CalledProcessError as e:
        print(f"Failed: {og_id}")
        with open(error_log, 'a') as logfile:
            logfile.write(f"Failed run {og_id}: {e}\n")

print("\nAll processing finished.")
