from Bio import SeqIO
import os

data_folder = "/nfs/scistore20/bronsgrp/jsaleh/data/missing_structures"
fasta_path = "/nfs/scistore20/bronsgrp/jsaleh/data/missing_sequences.fasta"
dat_path = os.path.join(data_folder, "missing_sequences_runs.dat")

with open(dat_path, 'w') as f:
    for record in SeqIO.parse(fasta_path, "fasta"):
        protein_id = record.id.split()[0]
        output_dir = f"{data_folder}/{protein_id}"
        a3m_file = f"{output_dir}/{protein_id}.a3m"
        f.write(f"{output_dir} {a3m_file}\n")
