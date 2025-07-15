"""
This script performs direct decoding (MSA-based consensus) for both filtered misleading and filtered normal reads
after double clustering. It uses MAFFT for multiple sequence alignment and evaluates reconstruction accuracy.
"""

import os
import shutil
import subprocess
from Bio import AlignIO
import Levenshtein

num = 10  # number of clusters
Lindex = 4

# Clear the given folder
def clear_folder(folder_path):
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print(f"Failed to delete {file_path}. Reason: {e}")

# Perform majority voting over aligned sequences
def Majority_Voting(maffres):
    linestr = ''
    baseWeightDict = {'A': 1, 'C': 1, 'G': 1, 'T': 1, '-': 0.4}
    weightVector = []
    for colid in range(len(maffres[0])):
        voterCounter = {k: 0 for k in baseWeightDict}
        for line in maffres:
            base = line[colid].upper()
            voterCounter[base] += baseWeightDict[base]
        cand = sorted(voterCounter.items(), key=lambda x: x[1], reverse=True)
        if cand[0][0] != '-':
            linestr += cand[0][0]
            weightVector.append(cand[0][1])
    return linestr, weightVector

# Get indices of the lowest tklast weights
def findMintklastIndex(weight, tklast):
    stdata = sorted(enumerate(weight), key=lambda x: x[1])
    return [i[0] for i in stdata[:tklast]]

# Save sequences to FASTA format
def write_sequences_to_fasta(filename, sequences):
    with open(filename, 'w', encoding='utf-8') as file:
        for index, sequence in enumerate(sequences):
            file.write(f'>Sequence{index + 1}\n{sequence}\n')

# Run MAFFT and return consensus
def mafft_alignment_and_consensus(cluster_file, fasta_file, output_file, mafft_path):
    with open(cluster_file, 'r') as infile:
        sequences = [line.strip() for line in infile]
    write_sequences_to_fasta(fasta_file, sequences)

    subprocess.run([mafft_path, "--auto", "--inputorder", "--quiet", "--thread", "-1", fasta_file],
                   stdout=open(output_file, "w"))
    alignment = AlignIO.read(output_file, "fasta")
    consensus, weight = Majority_Voting(alignment)

    # Adjust consensus length to 100 bp
    lklen = 100
    if len(consensus) > lklen:
        tklast = len(consensus) - lklen
        deleteindex = findMintklastIndex(weight, tklast)
        consensus = ''.join([base for i, base in enumerate(consensus) if i not in deleteindex])
    elif len(consensus) < lklen:
        consensus += 'A' * (lklen - len(consensus))
    return consensus

# Convert DNA to binary
def atcg_to_binary(atcg_sequences):
    binary_seq = ""
    for seq in atcg_sequences:
        for char in seq.strip():
            binary_seq += '0' if char in ['A', 'C'] else '1'
    return binary_seq

# Convert binary to text
def binary_to_text(binary_text):
    return ''.join([chr(int(binary_text[i:i+8], 2)) for i in range(0, len(binary_text), 8) if len(binary_text[i:i+8]) == 8])

# Compare similarity and error metrics
from scipy.spatial.distance import hamming

def compare_text(recovery_result, origtext):
    with open(recovery_result, 'r', encoding='utf-8') as file1, open(origtext, 'r', encoding='utf-8') as file2:
        file1_content = file1.read()
        file2_content = file2.read()
    compare_length = len(file2_content)
    file1_content = file1_content[:compare_length]
    file2_content = file2_content[:compare_length]
    edit_distance = Levenshtein.distance(file1_content, file2_content)
    max_distance = max(len(file1_content), len(file2_content))
    similarity_score = 1 - edit_distance / max_distance if max_distance > 0 else 0
    print('Similarity Score:', f'{similarity_score:.4f}')
    print('Edit Distance:', edit_distance)

def calculate_cer(original_text, recovered_text):
    truncated = recovered_text[:len(original_text)]
    edit_distance = Levenshtein.distance(original_text, truncated)
    cer = edit_distance / len(original_text) if len(original_text) > 0 else 0
    print(f"Character Error Rate (CER): {cer:.4f}")

def calculate_and_report_ber(file1, file2):
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        sequence1 = ''.join([line.strip() for line in f1])
        sequence2 = ''.join([line.strip() for line in f2])
    min_len = min(len(sequence1), len(sequence2))
    hamming_distance = hamming(list(sequence1[:min_len]), list(sequence2[:min_len])) * min_len
    ber = hamming_distance / min_len if min_len > 0 else 0
    print(f"Total bits: {min_len}")
    print(f"Error bits: {hamming_distance}")
    print(f"BER: {ber:.4f}")

# ---------------- Main process ----------------

def main():
    input_base_folder = "...\\decoding"
    subfolders = ["cluster_double_result_0", "cluster_double_result_1"]
    suffixes = ["0", "1"]
    mafft_path = "...\\decoding\\mafft-win\\mafft.bat"
    num = 10

    for subfolder, suffix in zip(subfolders, suffixes):
        folder_path = os.path.join(input_base_folder, subfolder)
        clear_folder("...\\decoding\\MSAcluster")
        consensus_sequences = []

        for k in range(num):
            cluster_file = os.path.join(folder_path, f'cluster_{k}_{suffix}.txt')
            fasta_file = f"...\\decoding\\MSAcluster\\MSAcluster_{k}_{suffix}.fasta"
            output_file = f"...\\decoding\\MSAcluster\\MSAcluster_{k}_{suffix}_result.txt"
            consensus = mafft_alignment_and_consensus(cluster_file, fasta_file, output_file, mafft_path)
            consensus_sequences.append(consensus)

            out_path = f"...\\decoding\\MSAcluster\\MSAcluster_{k}_{suffix}_consensus.txt"
            with open(out_path, 'w') as result_file:
                result_file.write(consensus)

        print(f"Consensus sequences generated for folder: {folder_path}")

        binary_sequences_with_index = []
        for k in range(num):
            ofilename = f"...\\decoding\\MSAcluster\\MSAcluster_{k}_{suffix}_consensus.txt"
            if os.path.exists(ofilename):
                with open(ofilename, "r") as ifile:
                    atcg_sequence = ifile.readline().strip()
                    print(f"ATCG sequence from file {k}: {atcg_sequence}")
                    binary_sequence = atcg_to_binary([atcg_sequence])
                    index = binary_sequence[:4]
                    sequence_content = binary_sequence[4:]
                    full_sequence = index + sequence_content
                    binary_sequences_with_index.append((index, full_sequence))
            else:
                print(f"File {ofilename} does not exist.")

        binary_sequences_with_index.sort(key=lambda x: x[0])

        with open('...\\decoding\\MSA_sorted_binary_sequences.txt', 'w') as ofile:
            for _, binary_seq in binary_sequences_with_index:
                ofile.write(binary_seq + '\n')

        binary_text = ''.join(seq[4:] for _, seq in binary_sequences_with_index)
        text_recovered = binary_to_text(binary_text)

        with open("...\\decoding\\MSA_recovered_message.txt", "w", encoding="utf-8") as output_file:
            output_file.write(text_recovered)

        # Evaluation
        file1 = '...\\decoding\\MSA_sorted_binary_sequences.txt'
        file2 = '...\\encoding\\original-01-message.txt'
        calculate_and_report_ber(file1, file2)

        original_text = open(r"...\\encoding\\original message.txt", 'r', encoding='utf-8').read()
        recovered_text = open(r"...\\decoding\\MSA_recovered_message.txt", 'r', encoding='utf-8').read()
        calculate_cer(original_text, recovered_text)
        compare_text(r"...\\decoding\\MSA_recovered_message.txt", r"...\\encoding\\original message.txt")

if __name__ == "__main__":
    main()
