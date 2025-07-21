import pandas as pd
from collections import defaultdict, Counter
import argparse

def hamming_distance(str1, str2):
    """
    Compute the Hamming distance between two strings.
    """
    if len(str1) != len(str2):
        raise ValueError("Strings must be of the same length to compute Hamming distance.")
    return sum(ch1 != ch2 for ch1, ch2 in zip(str1, str2))

def cluster_umis_with_freq_ref(dataframe, mismatch_threshold):
    """
    Clusters UMIs allowing for a certain number of mismatches.

    Args:
        dataframe (pd.DataFrame): Input dataframe with a column corresponding to UMI sequences.
        mismatch_threshold (int): Maximum allowed mismatches for UMIs to be grouped.

    Returns:
        pd.DataFrame: Original dataframe with an added column 'UMI_Group'.
    """
    # Get the list of UMIs
    umi_list = dataframe['Sequence'].tolist()
    umi_counts = Counter(umi_list) # count frequency for each UMI

    # Sort UMI by decreasing frequency
    sorted_umis = sorted(umi_counts.keys(), key=lambda umi: -umi_counts[umi])

    # Initialize variables for grouping
    umi_groups = defaultdict(list)  # To store grouped UMIs
    assigned_groups = {}  # To map each UMI to its group
    group_id = 0  # Group ID counter

    for umi in sorted_umis:
        if umi not in assigned_groups:
            # Create a new group for the current UMI
            group_id += 1
            umi_groups[group_id].append(umi)
            assigned_groups[umi] = group_id

            # Compare with all other UMIs to find matches within the threshold
            for other_umi in sorted_umis:
                if other_umi not in assigned_groups:
                    if hamming_distance(umi, other_umi) <= mismatch_threshold:
                        umi_groups[group_id].append(other_umi)
                        assigned_groups[other_umi] = group_id

    # Map the UMI to its group in the dataframe
    dataframe['UMI_group'] = dataframe.apply(lambda row: f"{assigned_groups[row['Sequence']]}_{row['seqnames.genome']}", axis=1)

    return dataframe

def main():
    parser = argparse.ArgumentParser(description="Cluster UMIs with a given mismatch threshold.")
    parser.add_argument("input_file", help="Path to the input file.")
    parser.add_argument("output_file", help="Path to the output file.")
    parser.add_argument("--mismatch_threshold", type=int, default=1, help="Maximum allowed mismatches for UMI grouping (default: 1).")

    args = parser.parse_args()

    # Read the input CSV file
    df = pd.read_csv(args.input_file, sep="\t")

    # Ensure the necessary columns are present
    if 'shearSite.genome' not in df.columns or 'Sequence' not in df.columns:
        raise ValueError("Input file must contain 'shearSite.genome' and 'Sequence' columns.")

    # Cluster UMIs
    df_clustered = cluster_umis_with_freq_ref(df, mismatch_threshold=args.mismatch_threshold)

    # Write the output to a CSV file
    df_clustered.to_csv(args.output_file, sep="\t", index=False)

if __name__ == "__main__":
    main()


