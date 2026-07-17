import os
import subprocess
import time

prot_dir = '/home/insan/python/new/dedup_proteins/final/modified/dNdS/Cleaned_Protein_FASTAs'
cds_dir  = '/home/insan/python/new/dedup_proteins/final/modified/dNdS/Cleaned_CDS_FASTAs'
out_dir  = '/home/insan/python/new/dedup_proteins/final/modified/dNdS/Codon_Alignments'
log_file = '/home/insan/python/new/dedup_proteins/final/modified/dNdS/pal2nal_errors.log'
threads  = 80

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

files = [f for f in os.listdir(prot_dir) if f.endswith('.fasta')]
files.sort()

with open(log_file, 'w') as log:
    for f in files:
        og_id = f.replace('.fasta', '')
        prot_in = os.path.join(prot_dir, f)
        cds_in  = os.path.join(cds_dir, f"{og_id}.cds.fasta")

        temp_aln = os.path.join(out_dir, f"{og_id}.prot.aln")
        final_nuc = os.path.join(out_dir, f"{og_id}.codon.nuc")

        if not os.path.exists(cds_in):
            continue

        subprocess.run(f"mafft --auto --thread {threads} {prot_in} > {temp_aln}", shell=True, check=True)
        time.sleep(0.1)

        # PAL2NAL
        cmd = f"pal2nal.pl {temp_aln} {cds_in} -output fasta -gc 1 -nomismatch -nogap > {final_nuc}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if os.path.exists(final_nuc) and os.path.getsize(final_nuc) == 0:
            error_msg = f"FAILED: {og_id} | Error: {result.stderr.strip()}\n"
            log.write(error_msg)
            print(error_msg.strip())
        else:

            if os.path.exists(temp_aln):
                os.remove(temp_aln)

print(f"\nDone. Check {log_file} for any remaining failures.")
