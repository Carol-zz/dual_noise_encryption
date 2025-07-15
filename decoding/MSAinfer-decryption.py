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

import Levenshtein
import numpy as np
import pickle
import subprocess
import os
import shutil
from Bio import AlignIO

num = 10  # number of clusters
Lindex = 4
classnum = num

def clear_directory(folder_path):
    """
    Clear all files and subfolders inside a given directory.
    Args:
        folder_path (str): Path to the directory to be cleared.
    """
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print(f'Failed to delete {file_path}. Reason: {e}')

# Convert DNA to carrier binary string: A/T -> 0, C/G -> 1
def atcg_to_carrier(basestring):
    return basestring.replace('A', '0').replace('a', '0').replace('T', '0').replace('t', '0') \
                     .replace('C', '1').replace('c', '1').replace('G', '1').replace('g', '1')

# Convert DNA (ACGT) to binary string: A/C -> 0, T/G -> 1
def atcg_to_binary(atcg_sequences):
    binary_seq = ""
    for seq in atcg_sequences:
        for char in seq.strip():
            if char in ['A', 'C']:
                binary_seq += '0'
            elif char in ['T', 'G']:
                binary_seq += '1'
    return binary_seq

# Convert binary string to ASCII text
def binary_to_text(binary_text):
    text = ""
    for i in range(0, len(binary_text), 8):
        byte = binary_text[i:i+8]
        if len(byte) == 8:
            text += chr(int(byte, 2))
    return text

# Find indices of the smallest `tklast` values in the weight vector
def findMintklastIndex(weight, tklast):
    stdata = sorted(enumerate(weight), key=lambda x: x[1])
    indexlist = [i[0] for i in stdata[:tklast]]
    return indexlist

# Perform majority voting over aligned sequences to generate a consensus
def Majority_Voting(maffres):
    linestr = ''
    baseWeightDict = {'A': 1, 'C': 1, '-': 0.4}
    voterCounter = {'A': 0, 'C': 0, '-': 0}
    weightVector = []
    for colid in range(len(maffres[0])):
        voterCounter = {'A': 0, 'C': 0, '-': 0}
        for line in maffres:
            base = line[colid].upper()
            if base in voterCounter:
                voterCounter[base] += baseWeightDict[base]
        cand = sorted(voterCounter.items(), key=lambda vt: vt[1], reverse=True)
        if cand[0][0] != '-':
            linestr += cand[0][0]
            weightVector.append(cand[0][1])
    return linestr, weightVector

# Use MAFFT to align sequences and generate consensus string of length `lklen`
def mafft_alignment_and_consensus(cluster_file, lklen=100):
    base_path = cluster_file.rsplit('.', 1)[0]
    fasta_file = f'{base_path}.fasta'
    output_file = f'{base_path}_consensus.fasta'

    with open(cluster_file, 'r') as infile, open(fasta_file, 'w') as outfile:
        for i, line in enumerate(infile):
            outfile.write(f'>Sequence{i + 1}\n{line.strip()}\n')

    mafft_path = ".../decoding/mafft-win/mafft.bat"
    subprocess.run([mafft_path, "--auto", "--inputorder", "--quiet", "--thread", "-1", fasta_file],
                   stdout=open(output_file, "w"))

    alignment = AlignIO.read(output_file, "fasta")
    aligned_seq, weight_vector = Majority_Voting(alignment)

    if len(aligned_seq) > lklen:
        delete_indices = findMintklastIndex(weight_vector, len(aligned_seq) - lklen)
        aligned_seq = ''.join([aligned_seq[i] for i in range(len(aligned_seq)) if i not in delete_indices])
    elif len(aligned_seq) < lklen:
        aligned_seq += 'A' * (lklen - len(aligned_seq))

    return aligned_seq

# Replace all 'A' with '1', and all 'C' with '0', and save to output file
def replace_ac_with_binary(input_path, output_path):
    try:
        with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
            for line in infile:
                modified_line = line.replace('A', '1').replace('C', '0')
                outfile.write(modified_line)
    except Exception as e:
        print(f"Error while processing file: {e}")


genecandstrlist = []

folder_path = ".../decoding/inferkey_cluster"
clear_directory(folder_path)

for i in range(classnum):
    input_filename = f".../decoding/clusterresult/cluster_{i}.txt"
    output_filename = f".../decoding/inferkey_cluster/inferkeycluster_{i}.txt"

    # Read, transform, and write sequences
    try:
        with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
            for line in infile:
                # First replace: A, T → 0; C, G → 1
                temp_line = line.replace('A', '0').replace('T', '0').replace('C', '1').replace('G', '1')
                # Second replace: 0 → C; 1 → A
                modified_line = temp_line.replace('1', 'A').replace('0', 'C')
                outfile.write(modified_line)
    except Exception as e:
        print(f"Error while processing file: {e}")

    # Perform MAFFT alignment and generate consensus key
    consensuskey = mafft_alignment_and_consensus(output_filename)
    genecandstrlist.append(consensuskey)

# Write inferred AC-style modukeys to file
output_path = ".../decoding/infered-modu-ACstr.txt"
with open(output_path, 'w') as outfile:
    for sequence in genecandstrlist:
        outfile.write(sequence + '\n')

# Generate final inferred modukey by second-round consensus
cluster_file = output_path
infermodukeyMAFFT = mafft_alignment_and_consensus(cluster_file)

# Convert final modukey from AC to binary
infermodukey = infermodukeyMAFFT.replace('A', '1').replace('C', '0')
print('inferred modukey is', infermodukey)

##########################################################################################
# Compare inferred modukey with the correct modukey (edit distance and hamming distance)
from scipy.spatial.distance import hamming

def read_key_file(filepath):
    key = []
    with open(filepath, "r") as file:
        for line in file:
            if line.strip():
                key.append(line.strip())
    if not key:
        raise ValueError("The key list is empty.")
    return key

key = read_key_file(".../key_generation/generated_key_2_5.txt")
print('the correct modukey is', key[0])

# Compare full inferred key with reference
inferkeyD = Levenshtein.distance(infermodukey[:len(key[0])], key[0])
print("Whole inferred Edit distance is:", inferkeyD)

inferkeyDH = hamming(list(infermodukey[:len(key[0])]), list(key[0])) * len(key[0])
print("Whole inferred Hamming distance is:", int(inferkeyDH))

##########################################################################################
# Compare each per-cluster inferred modukey to the ground truth key
def calculate_edit_distances(modu_file, key):
    with open(modu_file, "r") as file:
        line_number = 1
        for line in file:
            sequence = line.strip()
            if sequence:
                distance = Levenshtein.distance(sequence, key[0])
                print(f'Line {line_number}: Edit distance = {distance}')
                inferkeyH = hamming(list(sequence), list(key[0])) * len(key[0])
                print("  Inferred Hamming distance:", int(inferkeyH))
                print(sequence)
            line_number += 1

input_file = ".../decoding/infered-modu-ACstr.txt"
output_file = ".../decoding/infered-modu-01str.txt"
replace_ac_with_binary(input_file, output_file)
calculate_edit_distances(output_file, key)

# Store inferred modukey for later decoding
Modukey = infermodukey

# Copy original cluster sequences into filtered output directory
input_directory = ".../decoding/clusterresult"
output_directory = ".../decoding/inferkey_cluster"
os.makedirs(output_directory, exist_ok=True)

for filename in os.listdir(input_directory):
    if filename.endswith('.txt'):
        input_filepath = os.path.join(input_directory, filename)
        output_filepath = os.path.join(output_directory, f'inferkey_{filename}')
        with open(input_filepath, 'r') as file:
            content = file.read()
        with open(output_filepath, 'w') as f:
            f.write(content)

print("All sequences copied to inferkey_cluster directory.")

#######################################################################################
# Modulation-based error correction functions

def theta(a, b):
    if a == '-' or b == '-' or a != b:  # gap or mismatch
        return -1
    elif a == b:  # match
        return 1

def make_score_matrix(seq1, seq2):
    seq1 = '-' + seq1
    seq2 = '-' + seq2
    score_mat = {}
    trace_mat = {}

    for i, p in enumerate(seq1):
        score_mat[i] = {}
        trace_mat[i] = {}
        for j, q in enumerate(seq2):
            if i == 0:  # first row, gap in seq1
                score_mat[i][j] = -j
                trace_mat[i][j] = 1
                continue
            if j == 0:  # first column, gap in seq2
                score_mat[i][j] = -i
                trace_mat[i][j] = 2
                continue
            ul = score_mat[i - 1][j - 1] + theta(p, q)  # from up-left
            l = score_mat[i][j - 1] + theta('-', q)     # from left (gap in seq1)
            u = score_mat[i - 1][j] + theta(p, '-')     # from up (gap in seq2)
            picked = max([ul, l, u])
            score_mat[i][j] = picked
            trace_mat[i][j] = [ul, l, u].index(picked)  # record direction: 0=diag, 1=left, 2=up
    return score_mat, trace_mat

def traceback(seq1, seq2, trace_mat):
    seq1, seq2 = '-' + seq1, '-' + seq2
    i, j = len(seq1) - 1, len(seq2) - 1
    path_code = ''
    while i > 0 or j > 0:
        direction = trace_mat[i][j]
        if direction == 0:  # from diagonal (match or mismatch)
            i -= 1
            j -= 1
            path_code = '0' + path_code
        elif direction == 1:  # from left (gap in seq1)
            j -= 1
            path_code = '1' + path_code
        elif direction == 2:  # from up (gap in seq2)
            i -= 1
            path_code = '2' + path_code
    return path_code

def pretty_print_align(seq1, seq2, path_code):
    align1 = ''
    middle = ''
    align2 = ''
    for p in path_code:
        if p == '0':
            align1 += seq1[0]
            align2 += seq2[0]
            middle += '|' if seq1[0] == seq2[0] else ' '
            seq1 = seq1[1:]
            seq2 = seq2[1:]
        elif p == '1':
            align1 += '-'
            align2 += seq2[0]
            middle += ' '
            seq2 = seq2[1:]
        elif p == '2':
            align1 += seq1[0]
            align2 += '-'
            middle += ' '
            seq1 = seq1[1:]
    return (align1, align2)

def AlignTwoSeqs(seq1, seq2):
    score_mat, trace_mat = make_score_matrix(seq1, seq2)
    path_code = traceback(seq1, seq2, trace_mat)
    NewSeed, NewCandidate = pretty_print_align(seq1, seq2, path_code)
    return (NewSeed, NewCandidate)

def GroupReads(iflilename, multiplex):
    with open(iflilename, 'r', encoding='utf-8') as filept:
        data = filept.readlines()
        data = [vt[:-1] for vt in data]
        encode_num = len(data) // multiplex
        if len(data) % multiplex != 0:
            print("Simulation input file size error!")
        cls = [[] for _ in range(encode_num)]
        vcount = 0
        for eleid in range(encode_num):
            for tid in range(multiplex):
                cls[eleid].append(vcount)
                vcount += 1
        return cls, data

def ins_del_correct_v2(linestr, modustr):
    templatestr = modustr
    templineRealstrb = atcg_to_carrier(linestr)
    linestrb = atcg_to_binary(linestr)
    errorpos = []
    # align carrier bits with modulation key
    seq1, seq2 = AlignTwoSeqs(templineRealstrb, templatestr)
    resline = ""
    idvalue = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i] and seq1[i] != '-':
            resline += linestr[idvalue]
            errorpos.append('n')  # 'n' means normal match
            idvalue += 1
        elif seq1[i] == seq2[i] and seq1[i] == '-':
            continue
        elif seq1[i] != seq2[i] and seq1[i] == '-':  # deletion error
            errorpos.append('d')
            if seq2[i] == '0':
                resline += 'A' if np.random.rand() <= 0.5 else 'T'
            else:
                resline += 'G' if np.random.rand() <= 0.5 else 'C'
        elif seq1[i] != seq2[i] and seq2[i] == '-':  # insertion error
            idvalue += 1
        elif seq1[i] != seq2[i]:  # substitution error
            errorpos.append('s')
            if linestr[idvalue] in ['A', 'T']:
                resline += 'G' if np.random.rand() <= 0.5 else 'C'
            else:
                resline += 'A' if np.random.rand() <= 0.5 else 'T'
            idvalue += 1
    return resline, errorpos

def generateCandidate_1(newcontent, errorarray):
    tsize = len(newcontent)
    presline = ""
    for i in range(len(newcontent[0])):
        aplabet = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        for vx in range(tsize):
            if errorarray[vx][i] == 'n':
                aplabet[newcontent[vx][i]] += 1
        # pick the most frequent valid base
        presline += sorted(aplabet.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)[0][0]
    return presline

# Error correction of reads for each cluster and output consensus sequences
def cluster_indel_correction(clusterdict, data, modustr):
    rectiyreads = []
    errorstatarray = []
    # Perform indel correction for each read
    for dtline in data:
        tuele = ins_del_correct_v2(dtline, modustr)
        rectiyreads.append(tuele[0])
        errorstatarray.append(tuele[1])
    # Compute consensus sequence for each cluster
    content = []
    for clusele in clusterdict:
        if len(clusele) != 0:
            tmpcontent = [rectiyreads[ikd] for ikd in clusele]
            tmperray = [errorstatarray[ikd] for ikd in clusele]
            cite = generateCandidate_1(tmpcontent, tmperray)
            content.append(cite)
    return content

# Cluster reads and perform error correction
def cluster_majority(ifilename, ofilename, multiplex, modustr):
    cls, data = GroupReads(ifilename, multiplex)
    rectifyreads = cluster_indel_correction(cls, data, modustr)
    with open(ofilename, "w", encoding='utf-8') as wfile:
        for line in rectifyreads:
            wfile.write(line)
            wfile.write("\n")

def comFileLikely1(filename1):
    """Compare the character distribution between decoded file and original reference"""
    chrdict = dict()
    filesize = 0
    with open("template_refer_lib.lib", "rb") as file1:
        chrdict = pickle.load(file1)
        filesize = pickle.load(file1)
    with open(filename1, "rb") as file:
        data = file.read()
        contentfile = "".join([chr(b) for b in data])
    count = 0
    for vch in contentfile:
        if vch in chrdict and chrdict[vch] > 0:
            chrdict[vch] -= 1
        else:
            count += 1  # character not in reference
    total_diff = count + sum(abs(chrdict[vch]) for vch in chrdict if chrdict[vch] > 0)
    return 1 - total_diff / filesize

#################################################################################################################
import Levenshtein
from scipy.spatial.distance import hamming

def compare_text(recovery_result, origtext):
    """Compare similarity and edit distance between recovered and original text."""
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
    """Calculate Character Error Rate (CER) between original and recovered text."""
    truncated_recovered_text = recovered_text[:len(original_text)]
    edit_distance = Levenshtein.distance(original_text, truncated_recovered_text)
    cer = edit_distance / len(original_text) if len(original_text) > 0 else 0
    print(f"Character Error Rate (CER): {cer:.4f}")

def calculate_and_report_ber(file1, file2):
    """Calculate and report Bit Error Rate (BER) between two binary sequence files."""
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        sequence1 = ''.join([line.strip() for line in f1])
        sequence2 = ''.join([line.strip() for line in f2])

    min_len = min(len(sequence1), len(sequence2))
    hamming_distance = hamming(list(sequence1[:min_len]), list(sequence2[:min_len])) * min_len

    ber = hamming_distance / min_len if min_len > 0 else 0
    print(f"Total bits: {min_len}")
    print(f"Error bits: {hamming_distance}")
    print(f"BER: {ber:.4f}")

#########################################################################################

for k in range(num):
    input_file = f'.../decoding/inferkey_cluster/inferkey_cluster_{k}.txt'
    with open(input_file, 'r', encoding='utf-8') as filept:
        data = [line.strip() for line in filept]
        depth = len(data)
    # Perform correction and write consensus result
    cluster_majority(input_file, input_file.replace('.txt', '_consensus.txt'), depth, Modukey)

# Parse consensus results and build binary sequences with index
binary_sequences_with_index = []
for k in range(num):
    consensus_file = f'.../decoding/inferkey_cluster/inferkey_cluster_{k}_consensus.txt'
    if os.path.exists(consensus_file):
        with open(consensus_file, "r") as ifile:
            atcg_sequence = ifile.readline().strip()
            binary_sequence = atcg_to_binary(atcg_sequence)
            index = binary_sequence[:4]  # First 4 bits as index
            sequence_content = binary_sequence[4:]
            full_sequence = index + sequence_content
            binary_sequences_with_index.append((index, full_sequence))
    else:
        print(f"File {consensus_file} does not exist.")

# Sort by index and save the sorted binary sequences
binary_sequences_with_index.sort(key=lambda x: x[0])
sorted_binary_file = '.../decoding/infersorted_binary_sequences.txt'
with open(sorted_binary_file, 'w') as ofile:
    for index, binary_seq in binary_sequences_with_index:
        ofile.write(binary_seq + '\n')

# Convert all content (excluding index) to binary text and decode
binary_text = "".join(seq[4:] for _, seq in binary_sequences_with_index)
text_recovered = binary_to_text(binary_text)
with open('.../decoding/inferkey_recovered_message.txt', 'w', encoding='utf-8') as file:
    file.write(text_recovered)

###################################################################
# Evaluate final results
calculate_and_report_ber('.../decoding/infersorted_binary_sequences.txt',
                         '.../encoding/original-01-message.txt')

original_text = open('.../encoding/original message.txt', 'r', encoding='utf-8').read()
recovered_text = open('.../decoding/inferkey_recovered_message.txt', 'r', encoding='utf-8').read()
calculate_cer(original_text, recovered_text)

compare_text('.../decoding/inferkey_recovered_message.txt',
             '.../encoding/original message.txt')

