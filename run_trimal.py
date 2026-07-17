import os
import subprocess

ali_files='/home/insan/python/new/dedup_proteins/final/modified/mafftSingleCopy'
trim='/home/insan/python/new/dedup_proteins/final/modified/mafftSingleCopy/trimal'

if not os.path.exists(trim):
    os.makedirs(trim)

files=os.listdir(ali_files)
print(f"\nFound: {len(files)} files\n")

for f in files:
    if f.endswith('.fasta'):
        in_file=os.path.join(ali_files, f)
        out_name=f.rsplit('.', 1)[0] + '.fasta'
        out_file=os.path.join(trim, out_name)

        cmd=f"trimal -in {in_file} -out {out_file} -automated1"
        print(f"\nRunning trimal on {f}...\n")

        subprocess.run(cmd, shell=True)

print("\nAll trimming complete\n")
