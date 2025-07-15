"""
This script simulates the injection of artificial errors into modulated DNA sequences,
and subsequently simulates additional random channel noise as observed in sequencing.
It generates synthetic noisy DNA strands that can be used for evaluating the robustness
of modulation-based encryption in DNA storage.

Steps:
1. Load the modulation key and plaintext message.
2. Convert the message to binary, attach 4-bit indices.
3. Modulate the message using the key to generate DNA strands.
4. Inject artificial errors (replacement, insertion, deletion).
5. Apply random channel errors to simulate sequencing noise.
6. Save the final DNA strands for downstream decoding tests.
"""

import numpy as np
import random
import matplotlib.pyplot as plt
import argparse

# Parameters
Lpayload = 96
Lindex = 4
n = 4  # Number of modukeys per logical message (e.g., n=4 => 3 misleading + 1 correct)
Length = Lpayload + Lindex

# Command line arguments
parser = argparse.ArgumentParser(description="Simulate injected and channel errors.")
parser.add_argument('--error_rate', type=float, default=0.2, help='Artificial error injection rate')
parser.add_argument('--errorrate', type=float, default=0.2, help='Channel error rate')
args = parser.parse_args()
error_rate = args.error_rate
errorrate = args.errorrate
depth = 10  # sequencing depth per strand

# === Load modulation keys from file ===
def read_sequences(filename: str, n: int) -> list:
    """Read the first n binary sequences from file."""
    modukey = []
    with open(filename, 'r') as file:
        for i, line in enumerate(file):
            if i >= n:
                break
            modukey.append(line.strip())
    return modukey

# Load modulation keys
filename = r".../encoding/generated_key_2_5.txt"
modukey = read_sequences(filename, n)
for idx, key in enumerate(modukey):
    print(f"Sequence {idx + 1}: {key}")

# === Convert original message to binary with padding ===
with open(r".../encoding/original message.txt", "r", encoding="utf-8") as input_file:
    original_text = input_file.read()
binary_text = ''.join(format(ord(char), '08b') for char in original_text)
padding_length = (Lpayload - len(binary_text) % Lpayload) % Lpayload
padding = ''.join(random.choice(['0', '1']) for _ in range(padding_length))
binary_text += padding

# Add 4-bit index to each Lpayload block
unmoduMessage = []
for index in range(0, len(binary_text), Lpayload):
    binary_index = format(len(unmoduMessage), '04b').zfill(4)
    sequence_with_index = binary_index + binary_text[index:index + Lpayload]
    unmoduMessage.append(sequence_with_index)

# Save binary message with index to file
with open(r".../encoding/original-01-message.txt", "w") as output_file:
    for sequence in unmoduMessage:
        output_file.write(sequence + "\n")

# === Inject artificial group-wise errors ===
def introduce_errors_with_focus(sequences, error_rate):
    """
    Inject synthetic errors (replacement, insertion, deletion) at consistent positions
    across all misleading sequences in each group.
    """
    num_sequences = len(sequences)
    seq_length = len(sequences[0])
    num_errors = int(seq_length * error_rate / 100)

    error_indices = random.sample(range(seq_length), num_errors)
    print(error_indices)

    mutated_sequences = [list(seq) for seq in sequences]
    error_details = {'replacement': [], 'deletion': [], 'insertion': []}

    for i in range(num_sequences):
        if i % n != 0:  # skip reference (non-misleading) reads
            for idx in error_indices:
                error_type = random.choice(['replacement', 'deletion', 'insertion'])
                error_details[error_type].append((i, idx))
                if error_type == 'replacement':
                    mutated_sequences[i][idx] = random.choice(['A', 'T', 'C', 'G'])
                elif error_type == 'deletion':
                    mutated_sequences[i][idx] = ''
                elif error_type == 'insertion':
                    mutated_sequences[i].insert(idx + 1, random.choice(['A', 'T', 'C', 'G']))

    for i in range(len(mutated_sequences)):
        mutated_sequences[i] = ''.join([x for x in mutated_sequences[i] if x])  # remove empty slots

    return mutated_sequences, error_details

# === Modulate binary message with modulation key and inject errors ===
def generate_basestrings(Modukey, unmoduMessage, Length, error_rate):
    """
    Encode binary message using modulation key to DNA bases,
    and apply group-based injected errors.
    """
    basestrings = []
    for unmodu in unmoduMessage:
        group = []
        for moduk in Modukey:
            basestring = ''
            for i in range(Length):
                basestring += {
                    ('0', '0'): 'A',
                    ('0', '1'): 'T',
                    ('1', '0'): 'C',
                    ('1', '1'): 'G'
                }[(moduk[i], unmodu[i])]
            group.append(basestring)
        basestrings.extend(group)

    return introduce_errors_with_focus(basestrings, error_rate)

# Run and save modulated noisy strands
basestrings_with_errors, error_indices = generate_basestrings(modukey, unmoduMessage, Length, error_rate * 100)
output_file_path = r".../encoding/errorstrands_in.txt"
with open(output_file_path, "w") as output_file:
    for basestring in basestrings_with_errors:
        output_file.write(basestring + "\n")

# === Simulate sequencing-like random channel errors ===
def simulate_base_errors(ifilename, errorrate, multipleX, output_file):
    """
    Simulate random channel noise: substitutions, insertions, and deletions.
    Each sequence is replicated multipleX times with independent noise.
    """
    base_set = ['A', 'G', 'T', 'C']
    with open(ifilename, "r", encoding="UTF-8") as infile:
        sequences = infile.readlines()

    with open(output_file, "w", encoding="UTF-8") as outfile:
        for sequence in sequences:
            sequence = sequence.strip()
            for _ in range(multipleX):
                mutated_sequence = list(sequence)
                for i in range(len(mutated_sequence)):
                    if np.random.rand() <= errorrate:
                        error_type = np.random.choice(['substitution', 'insertion', 'deletion'])
                        if error_type == 'substitution':
                            alt = base_set[:]
                            alt.remove(mutated_sequence[i])
                            mutated_sequence[i] = random.choice(alt)
                        elif error_type == 'insertion':
                            mutated_sequence.insert(i, random.choice(base_set))
                        elif error_type == 'deletion':
                            mutated_sequence[i] = 'D'
                mutated_sequence = [b for b in mutated_sequence if b != 'D']
                outfile.write(''.join(mutated_sequence) + '\n')

# Apply channel noise
input_file = r".../encoding/errorstrands_in.txt"
output_file = r".../encoding/errorstrands_output.txt"
simulate_base_errors(input_file, errorrate, depth, output_file)

# === Optional: Visualize error positions on sequences ===
def plot_sequences_with_errors(sequences, error_indices):
    """
    Visualize the error types and their positions for debugging.
    Color codes:
        red: replacement
        blue: deletion
        green: insertion
    """
    num_sequences = len(sequences)
    seq_length = max(len(seq) for seq in sequences)
    print(error_indices)

    fig, ax = plt.subplots(figsize=(seq_length / 5, num_sequences))
    ax.set_title('Sequences with Error Positions Highlighted')
    ax.set_xlabel('Position')
    ax.set_ylabel('Sequence Index')

    for i, seq in enumerate(sequences):
        for j, char in enumerate(seq):
            color = 'black'
            weight = 'normal'
            if (i, j) in error_indices['replacement']:
                color = 'red'
                weight = 'bold'
            elif (i, j) in error_indices['deletion']:
                color = 'blue'
                weight = 'bold'
                char = '-'  # show deletion
            elif (i, j) in error_indices['insertion']:
                color = 'green'
                weight = 'bold'
            ax.text(j, i, char, va='center', ha='center', color=color, fontsize=10, fontweight=weight)

    ax.set_xlim(-1, seq_length)
    ax.set_ylim(-1, num_sequences)
    ax.set_xticks(range(seq_length))
    ax.set_yticks(range(num_sequences))
    ax.grid(True, color='lightgray')
    plt.show()

# Optional call to visualize:
# plot_sequences_with_errors(basestrings_with_errors, error_indices)
