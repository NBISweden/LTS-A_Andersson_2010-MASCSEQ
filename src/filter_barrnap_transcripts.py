from Bio import SeqIO
import gzip as gz
import pandas as pd


def main(sm):
    df = pd.read_csv(sm.input.barrnap, comment="#", sep="\t",
                     header=None, usecols=[0, 5, 8], index_col=0,
                     names=["transcript", "evalue", "name"])
    with gz.open(sm.output[0], 'wt') as fhout:
        for record in SeqIO.parse(gz.open(sm.input.ref, 'rt'), "fasta"):
            if record.id in df.index:
                continue
            fhout.write(f">{record.description}\n{record.seq}\n")


if __name__ == "__main__":
    main(snakemake)