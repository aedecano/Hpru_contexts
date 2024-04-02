import argparse
from Bio import SeqIO

def extract_sequences_in_range(genbank_file, start_product, end_product):
    sequences = {}
    start_found = False
    end_found = False
    for record in SeqIO.parse(genbank_file, "genbank"):
        sequences[record.id] = {}
        for feature in record.features:
            if 'product' in feature.qualifiers:
                product = feature.qualifiers['product'][0]
                if product == start_product:
                    start_found = True
                if start_found and not end_found:
                    if product not in sequences[record.id]:
                        sequences[record.id][product] = []
                    sequences[record.id][product].append(feature.extract(record.seq))
                if product == end_product:
                    end_found = True
                    break
    return sequences

def write_fasta_output(sequences, output_file):
    with open(output_file, "w") as f:
        for record_id, product_sequences in sequences.items():
            f.write(">{}\n".format(record_id))  # Write a single header for all sequences
            for product_name, sequences_list in product_sequences.items():
                for i, sequence in enumerate(sequences_list, start=1):
                    f.write(str(sequence) + "\n")

def main():
    parser = argparse.ArgumentParser(description="Extract sequences between start and end products from a GenBank file")
    parser.add_argument("genbank_file", help="Input GenBank file")
    parser.add_argument("-s", "--start_product", help="Start product name", required=True)
    parser.add_argument("-e", "--end_product", help="End product name", required=True)
    args = parser.parse_args()

    sequences = extract_sequences_in_range(args.genbank_file, args.start_product, args.end_product)
    output_file = args.genbank_file.replace(".gbk", "_cut.fasta")
    write_fasta_output(sequences, output_file)

if __name__ == "__main__":
    main()
