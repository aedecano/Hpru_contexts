import argparse
from Bio import SeqIO

def extract_sequences_by_product(genbank_file, product_names):
    sequences = {}
    for record in SeqIO.parse(genbank_file, "genbank"):
        sequences[record.id] = {}
        for feature in record.features:
            if 'product' in feature.qualifiers:
                product = feature.qualifiers['product'][0]
                if product in product_names:
                    if product not in sequences[record.id]:
                        sequences[record.id][product] = []
                    sequences[record.id][product].append(feature.extract(record.seq))
    return sequences

def write_fasta_output(sequences, output_file):
    with open(output_file, "w") as f:
        for record_id, product_sequences in sequences.items():
            for product_name, sequences_list in product_sequences.items():
                for i, sequence in enumerate(sequences_list, start=1):
                    f.write(">{}_{}_{}\n".format(record_id, product_name, i))
                    f.write(str(sequence) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Extract sequences by product from a GenBank file")
    parser.add_argument("genbank_file", help="Input GenBank file")
    parser.add_argument("-p", "--products", nargs="+", help="List of product names to extract", required=True)
    args = parser.parse_args()

    sequences = extract_sequences_by_product(args.genbank_file, args.products)
    output_file = args.genbank_file.replace(".gbk", "_cut.fasta")
    write_fasta_output(sequences, output_file)

if __name__ == "__main__":
    main()
