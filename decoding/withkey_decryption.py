"""
This script performs legitimate decryption after clustering in a DNA-based data storage system.
It filters reads based on edit distance to two modulation keys (one true, one misleading), selects
those closer to the correct key, and performs modulation decoding with error correction. The recovered
data is evaluated for accuracy using metrics such as BER (bit error rate), CER (character error rate),
edit distance, and similarity score compared to the original message.

Workflow:
1. Reads are clustered.
2. Sequences closer to the true modulation key are retained.
3. Modulation decoding and indel correction are applied.
4. Recovered binary sequences are assembled and converted to text.
5. Accuracy is measured against the original data.
"""

import os
import Levenshtein
import numpy as np
import pickle

# Parameters
n = 4               # modulation depth group size
num = 10            # number of clusters
Lindex = 4          # length of index bits

# Load true and misleading modulation keys (used for classification)
with open("D:\\data\\cl\\mutiple_confusion\\key_generation\\generated_key_2_5.txt", "r") as file:
    Modukeytrue = file.readline().strip()
    Modukeysimi = file.readline().strip()

# Convert DNA to binary carrier representation
def atcg_to_carrier(basestring):
    return basestring.replace('A', '0').replace('T', '0').replace('C', '1').replace('G', '1')

# Convert DNA to binary payload (bit string)
def atcg_to_binary(atcg_sequences):
    binary_seq = ""
    for seq in atcg_sequences:
        for char in seq.strip():
            if char in ['A', 'C']:
                binary_seq += '0'
            elif char in ['T', 'G']:
                binary_seq += '1'
    return binary_seq

# Convert binary to ASCII text
def binary_to_text(binary_text):
    text = ""
    for i in range(0, len(binary_text), 8):
        byte = binary_text[i:i+8]
        if len(byte) == 8:
            text += chr(int(byte, 2))
    return text

# Directory paths
input_directory = '...\\decoding\\clusterresult'  # (choose based on context)
output_directory = '...\\decoding\\withkey_cluster'

# Read standard ground-truth sequence labels
standard_cluster1_file = os.path.join('...\\decoding', 'standard_cluster1.txt')
standard_cluster2_file = os.path.join('...\\decoding', 'standard_cluster2.txt')
with open(standard_cluster1_file, 'r') as f:
    standard_cluster1 = set(f.read().splitlines())
with open(standard_cluster2_file, 'r') as f:
    standard_cluster2 = set(f.read().splitlines())

os.makedirs(output_directory, exist_ok=True)

# Filter reads based on distance to modulation keys
total_count_in_standard_cluster1 = 0
total_count_in_standard_cluster2 = 0

for filename in os.listdir(input_directory):
    if filename.endswith('.txt'):
        input_filepath = os.path.join(input_directory, filename)
        output_filepath = os.path.join(output_directory, f'withkey_{filename}')

        count_in_standard_cluster1 = 0
        count_in_standard_cluster2 = 0
        desired_sequences = []

        with open(input_filepath, 'r') as file:
            for line in file:
                atcg_sequence = line.strip()
                carrier_seq = atcg_to_carrier(atcg_sequence)
                dist_true = Levenshtein.distance(carrier_seq, Modukeytrue)
                dist_simi = Levenshtein.distance(carrier_seq, Modukeysimi)
                if dist_true <= dist_simi:
                    desired_sequences.append(atcg_sequence)

        for seq in desired_sequences:
            if seq in standard_cluster1:
                count_in_standard_cluster1 += 1
            elif seq in standard_cluster2:
                count_in_standard_cluster2 += 1

        total_count_in_standard_cluster1 += count_in_standard_cluster1
        total_count_in_standard_cluster2 += count_in_standard_cluster2

        print(f"File: {filename} - Filtered sequences: {len(desired_sequences)}")
        print(f"  From cluster1: {count_in_standard_cluster1}")
        print(f"  From cluster2: {count_in_standard_cluster2}")

        with open(output_filepath, 'w') as f:
            for seq in desired_sequences:
                f.write(seq + '\n')

total_desired_sequences = total_count_in_standard_cluster1 + total_count_in_standard_cluster2
overall_purity = (total_count_in_standard_cluster1 / total_desired_sequences) * 100 if total_desired_sequences > 0 else 0
print(f"Total cluster1: {total_count_in_standard_cluster1}")
print(f"Total cluster2: {total_count_in_standard_cluster2}")
print(f"Overall purity of selected reads: {overall_purity:.2f}%")
print("Filtering complete and saved to output directory.")

# ---------- Below: Modulation decoding and consensus correction utilities ----------
# (Includes sequence alignment, majority voting, error correction, etc.)
# -- omitted here for brevity, but will be included in GitHub version --

# ---------- Main Decryption Pipeline ----------
def main():
    skipped_k_values = []
    for k in range(num):
        ofilename = f'...\\decoding\\withkey_cluster\\withkey_cluster_{k}.txt'
        with open(ofilename, 'r', encoding='utf-8') as filept:
            data = [vt.strip() for vt in filept.readlines()]
            depth = len(data)

        if depth != 0:
            cluster_majority(ofilename, ofilename.replace('.txt', '_consensus.txt'), depth, Modukeytrue)
        else:
            skipped_k_values.append(k)
            consensus_filename = ofilename.replace('.txt', '_consensus.txt')
            with open(consensus_filename, 'w', encoding='utf-8') as file:
                file.write("")
            print(f"Cleared consensus for empty file {consensus_filename}")

    # Assemble binary sequences
    binary_sequences_with_index = []
    for k in range(num):
        if k not in skipped_k_values:
            ofilename = f'...\\decoding\\withkey_cluster\\withkey_cluster_{k}_consensus.txt'
            if os.path.exists(ofilename):
                with open(ofilename, "r") as ifile:
                    atcg_sequence = ifile.readline().strip()
                    binary_sequence = atcg_to_binary(atcg_sequence)
                    index = binary_sequence[:4]
                    sequence_content = binary_sequence[4:]
                    full_sequence = index + sequence_content
                    binary_sequences_with_index.append((index, full_sequence))

    binary_sequences_with_index.sort(key=lambda x: x[0])

    with open('...\\decoding\\sorted_binary_sequences.txt', 'w') as ofile:
        for _, binary_seq in binary_sequences_with_index:
            ofile.write(binary_seq + '\n')

    binary_text = "".join(seq[4:] for _, seq in binary_sequences_with_index)
    text_recovered = binary_to_text(binary_text)
    with open(r"...\\decoding\\withkey_recovered_message.txt", 'w', encoding='utf-8') as file:
        file.write(text_recovered)

    # Evaluation
    file1 = '...\\decoding\\sorted_binary_sequences.txt'
    file2 = '...\\encoding\\original-01-message.txt'
    calculate_and_report_ber(file1, file2)

    original_text = open(r"...\\encoding\\original message.txt", 'r', encoding='utf-8').read()
    recovered_text = open(r"...\\decoding\\withkey_recovered_message.txt", 'r', encoding='utf-8').read()

    calculate_cer(original_text, recovered_text)
    compare_text(
        r"...\\decoding\\withkey_recovered_message.txt",
        r"...\\encoding\\original message.txt"
    )

if __name__ == "__main__":
    main()
