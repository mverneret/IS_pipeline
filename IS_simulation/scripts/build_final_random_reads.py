import sys
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def random_replace(template: str, replace_dict: dict) -> str:
    return ''.join(random.choice(replace_dict.get(c, c)) for c in template)

def main(flanking_fasta, linkers_fasta, output_fasta):
    # Load linker + LTR sequences into a dictionary
    linker_dict = SeqIO.to_dict(SeqIO.parse(linkers_fasta, "fasta"))

    # Prepare output list
    final_records = []

    for record in SeqIO.parse(flanking_fasta, "fasta"):
        name = record.id
        flank_seq = str(record.seq)
        direction = 'sens' if random.random() < 0.5 else 'antisens'
        is_5prime = '5prime' in name
        seq_parts = []

        if is_5prime:
            random_V = random_replace("TTTVVVVTTVVVVTTVVVVTTVVVVTTT", {'V': ['A', 'C', 'G']})
            seq_parts.append(str(linker_dict['Linker_sens_before'].seq))
            seq_parts.append(random_V)
            seq_parts.append(str(linker_dict['Linker_sens_after'].seq))
            seq_parts.append(flank_seq)

            if direction == 'sens':
                seq_parts.append(str(linker_dict['LTR5_sens'].seq))
            else:
                seq_parts.append(str(linker_dict['LTR3_RC'].seq))

        else:  # 3prime
            random_B = random_replace("AAABBBBAABBBBAABBBBAABBBBAAA", {'B': ['T', 'C', 'G']})

            if direction == 'sens':
                seq_parts.append(str(linker_dict['LTR3_sens'].seq))
            else:
                seq_parts.append(str(linker_dict['LTR5_RC'].seq))

            seq_parts.append(flank_seq)
            seq_parts.append(str(linker_dict['Linker_RC_before'].seq))
            seq_parts.append(random_B)
            seq_parts.append(str(linker_dict['Linker_RC_after'].seq))

        full_seq = ''.join(seq_parts)

        # Random strand orientation for output
        output_strand = random.choice(['+', '-'])
        final_seq = Seq(full_seq)
        if output_strand == '-':
            final_seq = final_seq.reverse_complement()

        # Add final sequence with proper name
        new_id = f"{name}_{direction}_{output_strand}"
        final_records.append(SeqRecord(final_seq, id=new_id, description=""))

    # Write output
    SeqIO.write(final_records, output_fasta, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 build_final_random_reads.py <flanking_fasta> <linkers_LTR_fasta> <output_fasta>")
        sys.exit(1)

    flanking_fasta = sys.argv[1]
    linkers_fasta = sys.argv[2]
    output_fasta = sys.argv[3]

    main(flanking_fasta, linkers_fasta, output_fasta)
