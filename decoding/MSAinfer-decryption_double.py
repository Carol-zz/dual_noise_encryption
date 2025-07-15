"""
This script implements the 'inferMSA_double' attack method:
It performs inference of the modulation key using MAFFT-based consensus sequences from clustered subsets
(e.g., 'filtered normal reads' and 'filtered misleading reads'), simulating a two-stage plaintext attack.
"""

import numpy as np
import pickle
import subprocess
from Bio import AlignIO
import os
import shutil
import Levenshtein
from scipy.spatial.distance import hamming

num = 10  # Number of clusters
Lindex = 4
classnum = num

def clear_directory(folder_path):
    """
    Clear all files and subfolders in the specified folder.
    Args:
        folder_path (str): path to the folder to clear.
    """
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print(f"Failed to delete {file_path}. Reason: {e}")

############################################################################

def atcg_to_carrier(basestring):
    """
    Convert DNA bases to observed carrier bits.
    A/T -> 0, C/G -> 1
    """
    return basestring.replace('A', '0').replace('a', '0')\
                     .replace('T', '0').replace('t', '0')\
                     .replace('C', '1').replace('c', '1')\
                     .replace('G', '1').replace('g', '1')

def atcg_to_binary(atcg_sequences):
    """
    Convert ATCG sequences into binary strings using A/C=0, T/G=1 rule.
    """
    binary_seq = ""
    for seq in atcg_sequences:
        for char in seq.strip():
            if char in ['A', 'C']:
                binary_seq += '0'
            elif char in ['T', 'G']:
                binary_seq += '1'
    return binary_seq

def binary_to_text(binary_text):
    """
    Convert binary string into ASCII text by 8-bit chunks.
    """
    text = ""
    for i in range(0, len(binary_text), 8):
        byte = binary_text[i:i+8]
        if len(byte) == 8:
            text += chr(int(byte, 2))
    return text

#################################################################################

def findMintklastIndex(weight, tklast):
    """
    Find indices of the smallest tklast weights.
    """
    stdata = sorted(enumerate(weight), key=lambda x: x[1])
    indexlist = [i[0] for i in stdata[:tklast]]
    return indexlist

def Majority_Voting(maffres):
    """
    Compute consensus sequence by weighted voting across aligned columns.
    Only bases A and C and gap '-' are considered.
    """
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

def mafft_alignment_and_consensus(cluster_file, lklen=100):
    """
    Perform multiple sequence alignment using MAFFT and extract the consensus sequence.
    The consensus is trimmed or padded to length lklen.
    """
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

def replace_ac_with_binary(input_path, output_path):
    """
    Convert 'A' to '1', 'C' to '0' in a file and write to another file.
    Args:
        input_path (str): input file path.
        output_path (str): output file path.
    """
    try:
        with open(input_path, 'r') as infile, open(output_path, 'w') as outfile:
            for line in infile:
                modified_line = line.replace('A', '1').replace('C', '0')
                outfile.write(modified_line)
    except Exception as e:
        print(f"Error while processing file: {e}")

##########################################################################################
# Compare inferred modulation key(s) with ground-truth key
def read_key_file(filepath):
    """
    Read ground-truth modulation key(s) from file.
    """
    key = []
    with open(filepath, "r") as file:
        for line in file:
            if line.strip():
                key.append(line.strip())
    if not key:
        raise ValueError("The key list is empty.")
    return key

def calculate_edit_distances(modu_file, key):
    """
    Compute edit and Hamming distances between each inferred key and the true key.
    """
    with open(modu_file, "r") as file:
        line_number = 1
        for line in file:
            sequence = line.strip()
            if sequence:
                distance = Levenshtein.distance(sequence, key[0])
                print(f'Line {line_number}: Edit distance = {distance}')
                inferkeyH = hamming(list(sequence), list(key[0])) * len(key[0])
                print("Inferred Hamming distance:", int(inferkeyH))
                print(sequence)
            line_number += 1

#######################################################################################
# Indel-aware modulation error correction based on alignment
def theta(a, b):
    """
    Scoring function:
    match = +1, mismatch or gap = -1
    """
    if a == '-' or b == '-' or a != b:
        return -1
    return 1

def make_score_matrix(seq1, seq2):
    """
    Construct dynamic programming score and trace matrices for alignment.
    """
    seq1 = '-' + seq1
    seq2 = '-' + seq2
    score_mat = {}
    trace_mat = {}

    for i, p in enumerate(seq1):
        score_mat[i] = {}
        trace_mat[i] = {}
        for j, q in enumerate(seq2):
            if i == 0:
                score_mat[i][j] = -j
                trace_mat[i][j] = 1  # gap in seq1
                continue
            if j == 0:
                score_mat[i][j] = -i
                trace_mat[i][j] = 2  # gap in seq2
                continue
            ul = score_mat[i - 1][j - 1] + theta(p, q)
            l = score_mat[i][j - 1] + theta('-', q)
            u = score_mat[i - 1][j] + theta(p, '-')
            picked = max([ul, l, u])
            score_mat[i][j] = picked
            trace_mat[i][j] = [ul, l, u].index(picked)  # direction: 0=ul, 1=l, 2=u
    return score_mat, trace_mat

def traceback(seq1, seq2, trace_mat):
    """
    Reconstruct alignment path from traceback matrix.
    Returns:
        path_code: string of 0 (diag), 1 (left), 2 (up)
    """
    seq1, seq2 = '-' + seq1, '-' + seq2
    i, j = len(seq1) - 1, len(seq2) - 1
    path_code = ''
    while i > 0 or j > 0:
        direction = trace_mat[i][j]
        if direction == 0:
            i -= 1
            j -= 1
            path_code = '0' + path_code
        elif direction == 1:
            j -= 1
            path_code = '1' + path_code
        elif direction == 2:
            i -= 1
            path_code = '2' + path_code
    return path_code

def pretty_print_align(seq1, seq2, path_code):
    """
    Build alignment display from path code.
    Returns:
        aligned_seq1, aligned_seq2
    """
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
    return align1, align2

def AlignTwoSeqs(seq1, seq2):
    """
    Align two binary strings and return aligned sequences.
    """
    score_mat, trace_mat = make_score_matrix(seq1, seq2)
    path_code = traceback(seq1, seq2, trace_mat)
    NewSeed, NewCandidate = pretty_print_align(seq1, seq2, path_code)
    return NewSeed, NewCandidate

def GroupReads(iflilename, multiplex):
    """
    Group reads into clusters for decoding.
    Args:
        iflilename: input file path.
        multiplex: number of reads per message.
    Returns:
        cls: list of cluster indices.
        data: list of reads.
    """
    with open(iflilename, 'r', encoding='utf-8') as filept:
        data = filept.readlines()
        data = [vt.strip() for vt in data]
        encode_num = len(data) // multiplex
        if len(data) % multiplex != 0:
            print("Simulation source file has error!")
        cls = [[] for _ in range(encode_num)]
        vcount = 0
        for eleid in range(encode_num):
            for tid in range(multiplex):
                cls[eleid].append(vcount)
                vcount += 1
        return cls, data


#################################################################################################################
# 准确率、错误率指标计算
def compare_text(recovery_result, origtext):
    """计算恢复文本与原始文本之间的相似度和编辑距离，并打印结果。"""
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
def ins_del_correct_v2(linestr, modustr):
    """
    Perform indel and substitution correction on a single read using modulation key.

    Args:
        linestr: Original DNA read.
        modustr: Inferred modulation key (binary string).

    Returns:
        resline: Corrected DNA read.
        errorpos: List of error types at each position (n=normal, d=deletion, s=substitution).
    """
    templatestr = modustr
    templineRealstrb = atcg_to_carrier(linestr)
    linestrb = atcg_to_binary(linestr)
    errorpos = []

    seq1, seq2 = AlignTwoSeqs(templineRealstrb, templatestr)
    resline = ""
    idvalue = 0

    for i in range(len(seq1)):
        if seq1[i] == seq2[i] and seq1[i] != '-':
            resline += linestr[idvalue]
            errorpos.append('n')  # normal match
            idvalue += 1
        elif seq1[i] == seq2[i] and seq1[i] == '-':
            continue
        elif seq1[i] != seq2[i] and seq1[i] == '-':  # deletion
            errorpos.append('d')
            if seq2[i] == '0':
                resline += 'A' if np.random.rand() <= 0.5 else 'T'
            else:
                resline += 'G' if np.random.rand() <= 0.5 else 'C'
        elif seq1[i] != seq2[i] and seq2[i] == '-':  # insertion
            idvalue += 1
        elif seq1[i] != seq2[i]:  # substitution
            errorpos.append('s')
            if linestr[idvalue] in ['A', 'T']:
                resline += 'G' if np.random.rand() <= 0.5 else 'C'
            else:
                resline += 'A' if np.random.rand() <= 0.5 else 'T'
            idvalue += 1

    return resline, errorpos

def generateCandidate_1(newcontent, errorarray):
    """
    Generate a consensus sequence based on corrected reads and error labels.

    Args:
        newcontent: List of corrected reads.
        errorarray: List of corresponding error labels.

    Returns:
        presline: Consensus DNA sequence.
    """
    tsize = len(newcontent)
    presline = ""

    for i in range(len(newcontent[0])):
        aplabet = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        for vx in range(tsize):
            if errorarray[vx][i] == 'n':
                aplabet[newcontent[vx][i]] += 1
        presline += sorted(aplabet.items(), key=lambda kv: (kv[1], kv[0]), reverse=True)[0][0]

    return presline

def cluster_indel_correction(clusterdict, data, modustr):
    """
    Perform indel correction and compute consensus sequences for each cluster.

    Args:
        clusterdict: List of clusters, each containing indices.
        data: List of reads.
        modustr: Inferred modulation key.

    Returns:
        content: List of consensus sequences for each cluster.
    """
    rectiyreads = []
    errorstatarray = []

    for dtline in data:
        corrected_read, error_labels = ins_del_correct_v2(dtline, modustr)
        rectiyreads.append(corrected_read)
        errorstatarray.append(error_labels)

    content = []
    for clusele in clusterdict:
        if len(clusele) != 0:
            tmpcontent = [rectiyreads[ikd] for ikd in clusele]
            tmperray = [errorstatarray[ikd] for ikd in clusele]
            cite = generateCandidate_1(tmpcontent, tmperray)
            content.append(cite)

    return content

def cluster_majority(ifilename, ofilename, multiplex, modustr):
    """
    Perform clustering and error correction, then output consensus reads.

    Args:
        ifilename: Input file with reads.
        ofilename: Output file for consensus reads.
        multiplex: Number of reads per message.
        modustr: Inferred modulation key.
    """
    cls, data = GroupReads(ifilename, multiplex)
    rectifyreads = cluster_indel_correction(cls, data, modustr)
    with open(ofilename, "w", encoding='utf-8') as wfile:
        for line in rectifyreads:
            wfile.write(line + "\n")

def comFileLikely1(filename1):
    """
    Compare character frequency between a decoded file and reference character distribution.

    Args:
        filename1: Path to decoded file.

    Returns:
        likelihood: Similarity score between file and reference (1.0 = identical).
    """
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
            count += 1

    diff_sum = count + sum(abs(chrdict[vch]) for vch in chrdict if chrdict[vch] > 0)
    return 1 - diff_sum / filesize

    def calculate_cer(original_text, recovered_text):
        """
        Calculate and print Character Error Rate (CER) between original and recovered texts.

        Args:
            original_text: Ground truth string.
            recovered_text: Decoded string to be evaluated.
        """
        truncated_recovered_text = recovered_text[:len(original_text)]
        edit_distance = Levenshtein.distance(original_text, truncated_recovered_text)
        cer = edit_distance / len(original_text) if len(original_text) > 0 else 0
        print(f"Character Error Rate (CER): {cer:.4f}")

    def calculate_and_report_ber(file1, file2):
        """
        Calculate and print Bit Error Rate (BER) using Hamming distance between two binary files.

        Args:
            file1: Path to recovered binary file.
            file2: Path to original binary file.
        """
        with open(file1, 'r') as f1, open(file2, 'r') as f2:
            sequence1 = ''.join([line.strip() for line in f1])
            sequence2 = ''.join([line.strip() for line in f2])

        min_len = min(len(sequence1), len(sequence2))
        hamming_distance = hamming(list(sequence1[:min_len]), list(sequence2[:min_len])) * min_len

        ber = hamming_distance / min_len if min_len > 0 else 0
        print(f"Total bits: {min_len}")
        print(f"Error bits: {hamming_distance}")
        print(f"BER: {ber:.4f}")

    def main(input_directory, file_suffix):
        """
        Perform modulation key inference and decoding based on double-clustered reads.

        Args:
            input_directory: Path to folder containing cluster_*_{suffix}.txt files.
            file_suffix: Suffix used to identify normal or misleading reads.
        """
        genecandstrlist = []
        folder_path = ".../decoding/inferkey_cluster"
        clear_directory(folder_path)

        # Step 1: Preprocess input clusters to generate intermediate AC-only reads
        for i in range(classnum):
            input_filename = os.path.join(input_directory, f"cluster_{i}_{file_suffix}.txt")
            output_filename = f".../decoding/inferkey_cluster/inferkeycluster_{i}.txt"
            try:
                with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
                    for line in infile:
                        # First replace: A, T → 0; C, G → 1
                        temp_line = line.replace('A', '0').replace('T', '0').replace('C', '1').replace('G', '1')
                        # Second replace: 0 → C; 1 → A (AC-only representation)
                        modified_line = temp_line.replace('1', 'A').replace('0', 'C')
                        outfile.write(modified_line)
            except Exception as e:
                print(f"Error processing file: {e}")

            # Step 2: Generate consensus key from MSA of each cluster
            consensuskey = mafft_alignment_and_consensus(output_filename)
            genecandstrlist.append(consensuskey)

        # Step 3: Save inferred candidate keys from each cluster
        output_path = ".../decoding/infered-modu-ACstr.txt"
        with open(output_path, 'w') as outfile:
            for sequence in genecandstrlist:
                outfile.write(sequence + '\n')

        # Step 4: Align all cluster keys again to get global modulation key
        cluster_file = output_path
        infermodukeyMAFFT = mafft_alignment_and_consensus(cluster_file)
        infermodukey = infermodukeyMAFFT.replace('A', '1').replace('C', '0')
        print('Inferred modulation key is:', infermodukey)

    ##########################################################################################
    # Compare inferred modulation key with the ground truth key
    key = read_key_file(".../key_generation/generated_key_2_5.txt")
    print('The correct modulation key is:', key[0])

    # Compute edit distance and hamming distance
    inferkeyD = Levenshtein.distance(infermodukey[:len(key[0])], key[0])
    print("Whole inferred Edit distance is:", inferkeyD)
    inferkeyDH = hamming(list(infermodukey[:len(key[0])]), list(key[0])) * len(key[0])
    print("Whole inferred Hamming distance is:", int(inferkeyDH))

    ##########################################################################################
    # Evaluate edit distance between each inferred key and the true key
    input_file = ".../decoding/infered-modu-ACstr.txt"
    output_file = ".../decoding/infered-modu-01str.txt"
    replace_ac_with_binary(input_file, output_file)
    calculate_edit_distances(output_file, key)

    Modukey = infermodukey

    # Ensure output folder exists
    output_directory = ".../decoding/inferkey_cluster"
    os.makedirs(output_directory, exist_ok=True)

    # Copy all cluster files into the output folder (relabeling)
    for filename in os.listdir(input_directory):
        if filename.endswith('.txt'):
            input_filepath = os.path.join(input_directory, filename)
            output_filepath = os.path.join(output_directory, f'inferkey_{filename}')
            with open(input_filepath, 'r') as file:
                content = file.read()
            with open(output_filepath, 'w') as f:
                f.write(content)
    print("All cluster files have been copied to the correction folder.")

    ###################################################
    # Perform consensus generation for each corrected cluster
    for k in range(num):
        ofilename = f".../decoding/inferkey_cluster/inferkey_cluster_{k}_{file_suffix}.txt"
        with open(ofilename, 'r', encoding='utf-8') as filept:
            data = filept.readlines()
            data = [vt.strip() for vt in data]
            depth = len(data)
        cluster_majority(ofilename, ofilename.replace('.txt', '_consensus.txt'), depth, Modukey)

    # Decode all consensus sequences into binary form
    binary_sequences_with_index = []
    for k in range(num):
        ofilename = f".../decoding/inferkey_cluster/inferkey_cluster_{k}_{file_suffix}_consensus.txt"
        if os.path.exists(ofilename):
            with open(ofilename, "r") as ifile:
                atcg_sequence = ifile.readline().strip()
                binary_sequence = atcg_to_binary(atcg_sequence)
                index = binary_sequence[:4]
                sequence_content = binary_sequence[4:]
                full_sequence = index + sequence_content
                binary_sequences_with_index.append((index, full_sequence))
        else:
            print(f"File {ofilename} does not exist.")

    # Sort sequences by index and write them to file
    binary_sequences_with_index.sort(key=lambda x: x[0])
    output_filename = ".../decoding/infersorted_binary_sequences.txt"
    with open(output_filename, 'w') as ofile:
        for index, binary_seq in binary_sequences_with_index:
            ofile.write(binary_seq + '\n')

    # Convert binary string to text
    binary_text = "".join(seq[4:] for _, seq in binary_sequences_with_index)
    text_recovered = binary_to_text(binary_text)
    with open(".../decoding/inferkey_recovered_message.txt", 'w', encoding='utf-8') as file:
        file.write(text_recovered)

    ###################################################################
    # Evaluate BER and CER between original and recovered messages
    file1 = ".../decoding/infersorted_binary_sequences.txt"
    file2 = ".../encoding/original-01-message.txt"
    calculate_and_report_ber(file1, file2)

    original_text = open(".../encoding/original message.txt", 'r', encoding='utf-8').read()
    recovered_text = open(".../decoding/inferkey_recovered_message.txt", 'r', encoding='utf-8').read()
    calculate_cer(original_text, recovered_text)

    compare_text(".../decoding/inferkey_recovered_message.txt", ".../encoding/original message.txt")
