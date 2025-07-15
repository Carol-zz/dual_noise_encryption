"""
This script generates binary modulation keys that satisfy biochemical constraints,
such as limiting the maximum number of consecutive 0s or 1s.

It also generates misleading keys by mutating the original ones within a specified
Hamming distance range (min_dist ~ max_dist). One selected misleading key can be
duplicated multiple times to simulate multiple misleading reads.

Use case: DNA storage security experiments requiring confusion sequences.
"""

import random
from typing import List
import pandas as pd
from Crypto.Random import get_random_bytes

# === Generate random binary sequence with max run-length constraint ===
def generate_constrained_random_sequence(n: int, max_consecutive: int = 3) -> List[int]:
    """
    Generate a binary sequence of length n with at most max_consecutive identical bits.
    """
    sequence = []
    last_value = None
    consecutive_count = 0
    required_bytes = (n + 7) // 8
    random_bytes = get_random_bytes(required_bytes)

    for byte in random_bytes:
        for i in range(8):
            if len(sequence) == n:
                break
            value = (byte >> i) & 1
            if consecutive_count >= max_consecutive:
                value = 1 - last_value
                consecutive_count = 1
            elif last_value is None or value == last_value:
                consecutive_count += 1
            else:
                consecutive_count = 1
            sequence.append(value)
            last_value = value
    return sequence

# === Hamming distance computation ===
def hamming_distance(seq1: List[int], seq2: List[int]) -> int:
    return sum(el1 != el2 for el1, el2 in zip(seq1, seq2))

# === Mutate a binary sequence to ensure constrained distance ===
def mutate_constrained_sequence(base_seq: List[int], min_dist: int, max_dist: int, max_consecutive: int = 3) -> List[int]:
    """
    Generate a new sequence by flipping bits in base_seq, maintaining distance and run-length constraints.
    """
    n = len(base_seq)
    num_changes = random.randint(min_dist, max_dist)
    indices = random.sample(range(n), num_changes)
    new_seq = base_seq[:]

    for idx in indices:
        original_value = new_seq[idx]
        new_value = 1 - original_value
        if idx > 0 and idx < n - 1:
            if new_seq[idx - 1] == new_seq[idx + 1] == new_value:
                continue  # skip to avoid creating long runs
        new_seq[idx] = new_value
    return new_seq

# === Validate sequence distance from existing keys ===
def is_valid_sequence(new_seq: List[int], sequences: List[List[int]], min_dist: int, max_dist: int) -> bool:
    return all(min_dist <= hamming_distance(new_seq, seq) <= max_dist for seq in sequences)

# === Generate N modulation keys and misleading ones ===
def generate_sequences(num_sequences: int, sequence_length: int, min_dist: int, max_dist: int, max_consecutive: int = 3) -> List[List[int]]:
    sequences = [generate_constrained_random_sequence(sequence_length, max_consecutive)]
    while len(sequences) < num_sequences:
        new_seq = None
        attempts = 0
        while new_seq is None or not is_valid_sequence(new_seq, sequences, min_dist, max_dist):
            new_seq = mutate_constrained_sequence(sequences[-1], min_dist, max_dist, max_consecutive)
            attempts += 1
            if attempts > 100:
                new_seq = generate_constrained_random_sequence(sequence_length, max_consecutive)
                attempts = 0
        sequences.append(new_seq)
    return sequences

# === Display pairwise Hamming distance matrix ===
def display_distance_matrix(sequences: List[List[int]]) -> pd.DataFrame:
    num_sequences = len(sequences)
    distance_matrix = pd.DataFrame(index=range(1, num_sequences + 1), columns=range(1, num_sequences + 1))
    for i in range(num_sequences):
        for j in range(num_sequences):
            if i >= j:
                distance_matrix.at[i + 1, j + 1] = "-"
            else:
                dist = hamming_distance(sequences[i], sequences[j])
                distance_matrix.at[i + 1, j + 1] = dist
                distance_matrix.at[j + 1, i + 1] = dist
    return distance_matrix

# === Write sequences to file, with optional repeated entry ===
def write_sequences_to_file(sequences: List[List[int]], filename: str, repeat_index: int = None, repeat_times: int = 1):
    """
    Write sequences to file. Optionally duplicate one misleading key multiple times.

    Args:
        sequences: list of binary sequences
        filename: output file path
        repeat_index: index of misleading sequence to duplicate
        repeat_times: how many times to duplicate
    """
    with open(filename, 'w') as file:
        for i, seq in enumerate(sequences):
            if i == repeat_index:
                for _ in range(repeat_times):
                    file.write(''.join(map(str, seq)) + '\n')
            else:
                file.write(''.join(map(str, seq)) + '\n')

# === Example Usage ===
num_sequences = 2            # total number of keys (1 true + 1 misleading)
sequence_length = 100        # key length
min_dist = 20                # minimum edit distance between keys
max_dist = 20                # maximum edit distance between keys
max_consecutive = 2          # max allowed run-length (for biochemical constraint)

# Generate modulation key + misleading keys
sequences = generate_sequences(num_sequences, sequence_length, min_dist, max_dist, max_consecutive)

# Write keys to file, repeating the misleading one (e.g. 7 misleading reads)
write_sequences_to_file(sequences, "generated_key_2_20.txt", repeat_index=1, repeat_times=7)

# Print results
distance_matrix = display_distance_matrix(sequences)
for i, seq in enumerate(sequences):
    print(f"Sequence {i + 1}: {''.join(map(str, seq))}")
print("\nEdit Distance Matrix:")
print(distance_matrix)
