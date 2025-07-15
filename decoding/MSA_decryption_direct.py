"""
This script simulates a direct chosen-plaintext attack on DNA-based storage by performing
Multiple Sequence Alignment (MSA) using MAFFT over the entire pool of sequencing reads.
No clustering is used; instead, the global alignment result is used to generate a consensus
sequence via majority voting. The recovered binary sequence is then decoded to plaintext and
evaluated using standard similarity metrics (BER, CER, edit distance).
"""

import os
import shutil
import subprocess
from Bio import AlignIO
import Levenshtein

num = 10  # Number of clusters
Lindex = 4  # Index length

# Utility: Clear all files in a folder
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

# Perform majority voting from MSA result
def Majority_Voting(maffres):
    linestr = ''
    baseWeightDict = {'A': 1, 'C': 1, 'G': 1, 'T': 1, '-': 0.4}
    voterCounter = {'A': 0, 'C': 0, 'G': 0, 'T': 0, '-': 0}
    weigthVector = []
    for colid in range(len(maffres[0])):
        voterCounter = {key: 0 for key in baseWeightDict}
        for line in maffres:
            base = line[colid].upper()
            voterCounter[base] += baseWeightDict[base]
        cand = sorted(voterCounter.items(), key=lambda x: x[1], reverse=True)
        if cand[0][0] != '-':
            linestr += cand[0][0]
            weigthVector.append(cand[0][1])
    return linestr, weigthVector

# Return indexes of tklast smallest weights
def findMintklastIndex(weight, tklast):
    stdata = sorted(enumerate(weight), key=lambda x: x[1])
    indexlist = [i[0] for i in stdata]
    return indexlist[:tklast]

# Write a list of sequences to FASTA format
def write_sequences_to_fasta(filename, sequences):
    with open(filename, 'w', encoding='utf-8') as file:
        for index, sequence in enumerate(sequences):
            file.write(f'>Sequence{index + 1}\n{sequence}\n')

# Run MAFFT and extract consensus sequence
def mafft_alignment_and_consensus(cluster_file, fasta_file, output_file, mafft_path):
    with open(cluster_file, 'r') as infile:
        sequences = [line.strip() for line in infile]
    write_sequences_to_fasta(fasta_file, sequences)

    subprocess.run([mafft_path, "--auto", "--inputorder", "--quiet", "--thread", "-1", fasta_file],
                   stdout=open(output_file, "w"))

    alignment = AlignIO.read(output_file, "fasta")
    consensus, weight = Majority_Voting(alignment)

    # Adjust length to 100bp
    lklen = 100
    if len(consensus) > lklen:
        tklast = len(consensus) - lklen
        linestr = ""
        deleteindex = findMintklastIndex(weight, tklast)
        for vid, base in enumerate(consensus):
            if vid not in deleteindex:
                linestr += base
        consensus = linestr
    elif len(consensus) < lklen:
        consensus += 'A' * (lklen - len(consensus))
    return consensus

# ---------------- Binary and decoding utilities ------------------

def atcg_to_binary(atcg_sequences):
    binary_seq = ""
    for seq in atcg_sequences:
        for char in seq.strip():
            if char in ['A', 'C']:
                binary_seq += '0'
            elif char in ['T', 'G']:
                binary_seq += '1'
    return binary_seq

def binary_to_text(binary_text):
    text = ""
    for i in range(0, len(binary_text), 8):
        byte = binary_text[i:i+8]
        if len(byte) == 8:
            text += chr(int(byte, 2))
    return text

def recoveredDNA2text(ifilename, ofilename):
    with open(ifilename, "r") as ifile, open(ofilename, 'w') as ofile:
        atcg_sequences = [line.strip() for line in ifile]
        print("ATCG sequences:", atcg_sequences)
        binary_str = atcg_to_binary(atcg_sequences)
        print("Binary string:", binary_str)
        text = binary_to_text(binary_str)
        print("Recovered text:", text)
        ofile.write(text)

# ---------------- Evaluation Metrics ------------------

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
    truncated_recovered_text = recovered_text[:len(original_text)]
    edit_distance = Levenshtein.distance(original_text, truncated_recovered_text)
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

# ---------------- Main pipeline ------------------

def main():
    folder_path = "...\\decoding\\MSAcluster\\"
    mafft_path = "...\\decoding\\mafft-win\\mafft.bat"
    num = 10

    clear_folder(folder_path)
    consensus_sequences = []

    for k in range(num):
        cluster_file = f'...\\decoding\\clusterresult\\cluster_{k}.txt'
        fasta_file = f"...\\decoding\\MSAcluster\\MSAcluster_{k}.fasta"
        output_file = f"...\\decoding\\MSAcluster\\MSAcluster_{k}_result.txt"
        consensus = mafft_alignment_and_consensus(cluster_file, fasta_file, output_file, mafft_path)
        consensus_sequences.append(consensus)

        output_file_path = f"...\\decoding\\MSAcluster\\MSAcluster_{k}_consensus.txt"
        with open(output_file_path, 'w') as result_file:
            result_file.write(consensus)

    print('Consensus generation complete.')

    binary_sequences_with_index = []
    for k in range(num):
        ofilename = f'...\\decoding\\MSAcluster\\MSAcluster_{k}_consensus.txt'
        if os.path.exists(ofilename):
            with open(ofilename, "r") as ifile:
                atcg_sequence = ifile.readline().strip()
                print(f"ATCG sequence from file {k}:", atcg_sequence)
                binary_sequence = atcg_to_binary(atcg_sequence)
                index = binary_sequence[:4]
                sequence_content = binary_sequence[4:]
                full_sequence = index + sequence_content
                binary_sequences_with_index.append((index, full_sequence))
        else:
            print(f"File {ofilename} does not exist.")

    binary_sequences_with_index.sort(key=lambda x: x[0])

    output_filename = '...\\decoding\\MSA_sorted_binary_sequences.txt'
    with open(output_filename, 'w') as ofile:
        for index, binary_seq in binary_sequences_with_index:
            ofile.write(binary_seq + '\n')

    binary_text = "".join(seq[4:] for _, seq in binary_sequences_with_index)
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
