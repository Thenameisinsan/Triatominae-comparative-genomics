import sys,os,glob
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import ast
import re
import subprocess
import json


in_dir='/home/insan/python/old_proteins'
out_dir='/home/insan/python/new/dedup_proteins'
length_dir='/home/insan/python/new/lengths'
formatted='/home/insan/python/new/formatted'
deisoformed='/home/insan/python/new/de_isoformed'

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

if not os.path.exists(length_dir):
    os.makedirs(length_dir)

if not os.path.exists(formatted):
    os.makedirs(formatted)

if not os.path.exists(deisoformed):
    os.makedirs(deisoformed)

files = glob.glob(os.path.join(in_dir, "*.fasta"))
print(f"Found {len(files)} files processing...")

def create_file(out_):
    with open(out_, 'w') as out:
        out.write('')
        out.close()

for f in files:
    fname = os.path.basename(f)
    base_name = os.path.splitext(fname)[0]

    output = os.path.join(formatted, f"{base_name}_formatted.fasta")
    output_tsv = os.path.join(length_dir, f"{base_name}_lengths.tsv")
    output_3 = os.path.join(deisoformed, f"{base_name}_isoform_cleaned.fasta")
    output_4 = os.path.join(out_dir, f"{base_name}_dedup.fasta")

    create_file(output)
    create_file(output_tsv)
    create_file(output_3)
    create_file(output_4)

    print(f"processing {fname}...")

    for record in SeqIO.parse(f, 'fasta'):
        sequ=record.seq
        contig=record.id
        full_name=record.description
        trim=full_name.replace(contig, ' ').strip()
        newId=trim.split()
        newName='_'.join(newId)
        newRec=SeqRecord(id=newName, seq=sequ, description='')
        with open(output, 'a') as out1:
            SeqIO.write(newRec, out1, 'fasta')

        length=[newName, len(sequ)]
        with open(output_tsv, 'a') as out2:
            out2.write(f"{length[0]}\t{length[1]}\n")


    df=pd.read_csv(output_tsv, sep='\t', header=None)
    df.columns=['id', 'length']

    df['similar_names']=df['id'].str.split('_isoform').str[0]
    gru=df.groupby('similar_names')

    ls_id = []

    for group,info in gru:
        fixed=info.sort_values('length', ascending=False).drop_duplicates('similar_names', keep='first')
        ls_id.append(fixed.iloc[0,0])

    prot=0

    fasta2=SeqIO.parse(output, 'fasta')

    for i, rec in enumerate(fasta2):
        if rec.id in ls_id:
            new_recId=str(rec.id).replace("/", "_").replace(",", "_").replace("-", "_")
            new_rec=SeqRecord(id=new_recId, seq=rec.seq, description='')
            prot+=1
            with open(output_3, 'a') as out3:
                SeqIO.write(new_rec, out3, 'fasta')
    print(f"proteins written: {prot}")

    print(f"Running seqkit on {output_3}...")
    cmd = f"seqkit rmdup -s {output_3} > {output_4}"
    subprocess.run(cmd, shell=True, check=True)

print("All done")
