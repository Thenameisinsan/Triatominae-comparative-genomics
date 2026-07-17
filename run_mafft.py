fasta_files='/home/insan/python/new/dedup_proteins/final/modified/singleCopyOrthologs'
ali_files='/home/insan/python/new/dedup_proteins/final/modified/mafftSingleCopy'
threads=80

if not os.path.exists(ali_files):
    os.makedirs(ali_files)

files=os.listdir(fasta_files)
print(f'\nfound: {len(files)} files\n')

for f in files:
    if f.endswith('.fasta'):
        in_path=os.path.join(fasta_files, f)
        out_name=f.rsplit('.', 1)[0] + '.fasta'
        out_path=os.path.join(ali_files, out_name)

        cmd=f"mafft --auto --thread {threads} {in_path} > {out_path}"
        print(f"\nAligning {f}...\n")

        subprocess.run(cmd, shell=True)

print("\nAll alignments complete\n")
