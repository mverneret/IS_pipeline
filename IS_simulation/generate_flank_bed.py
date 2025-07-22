import sys
import random

def main(bed_in, bed_out, flank_len):
    with open(bed_in, 'r') as infile, open(bed_out, 'w') as outfile:
        for line in infile:
            if line.strip() == "":
                continue
            parts = line.strip().split()
            if len(parts) < 3:
                continue  # skip malformed lines

            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            insertion_id = f"{chrom}_{start}_{end}"

            direction = random.choice(["5prime", "3prime"])
            if direction == "5prime":
                flank_start = max(0, start - flank_len)
                flank_end = start
            else:
                flank_start = end
                flank_end = end + flank_len

            name = f"{insertion_id}_{direction}"
            score = 0
            strand = "+"

            outfile.write(f"{chrom}\t{flank_start}\t{flank_end}\t{name}\t{score}\t{strand}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python generate_flank_bed.py <input_bed> <output_bed> <flank_length>")
        sys.exit(1)

    input_bed = sys.argv[1]
    output_bed = sys.argv[2]
    flank_length = int(sys.argv[3])

    main(input_bed, output_bed, flank_length)
