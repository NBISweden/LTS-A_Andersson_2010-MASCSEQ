from Bio import SeqIO
import gzip as gz

def main(sm):
    rRNA = []
    for record in SeqIO.parse(sm.input.barrnap, "fasta"):
        rRNA.append(record.id)
    with gz.open(sm.output[0], 'wt') as fhout:
        for record in SeqIO.parse(gz.open(sm.input.ref, 'rt'), "fasta"):
            if record.id in rRNA:
                continue
            fhout.write(f">{record.description}\n{record.seq}\n")



if __name__ == "__main__":
    main(snakemake)