import argparse
import random
from Bio import SeqIO

def get_chromosome_lengths(fasta_file):
    chr_lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        chr_lengths[record.id] = len(record.seq)
    return chr_lengths

def pick_random_sites(chr_lengths, num_sites):
    chromosomes = list(chr_lengths.keys())
    lengths = list(chr_lengths.values())
    total_length = sum(lengths)
    weights = [length / total_length for length in lengths]

    sites = []
    for _ in range(num_sites):
        chosen_chr = random.choices(chromosomes, weights=weights, k=1)[0]
        max_pos = chr_lengths[chosen_chr]
        pos = random.randint(0, max_pos - 1)
        sites.append((chosen_chr, pos, pos + 1, '+'))
    return sites

def write_bed_file(sites, output_bed):
    with open(output_bed, 'w') as out:
        for chr, start, end, strand in sites:
            out.write(f"{chr}\t{start}\t{end}\t.\t0\t{strand}\n")

def main():
    parser = argparse.ArgumentParser(description="Generate random integration sites from a genome FASTA file.")
    parser.add_argument('--fasta', required=True, help='Input FASTA file of the genome assembly')
    parser.add_argument('--output', required=True, help='Output BED file with random integration sites')
    parser.add_argument('--num-sites', type=int, required=True, help='Number of random integration sites to generate')

    args = parser.parse_args()

    chr_lengths = get_chromosome_lengths(args.fasta)
    random_sites = pick_random_sites(chr_lengths, args.num_sites)
    write_bed_file(random_sites, args.output)

if __name__ == "__main__":
    main()
