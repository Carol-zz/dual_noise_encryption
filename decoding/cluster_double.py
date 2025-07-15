"""
This script implements the first step of decryption for a DNA-based storage encryption system.
It clusters the sequencing reads into 10 groups (each corresponding to an original message block),
and performs a second-level binary clustering (2 clusters) on each group to simulate
ciphertext-only attacks using various clustering algorithms (e.g., K-means, hierarchical, spectral).
"""

import os
os.environ['OMP_NUM_THREADS'] = '1'
import matplotlib.pyplot as plt
import numpy as np
import Levenshtein
from sklearn.cluster import AgglomerativeClustering
import networkx as nx

depth = 10
n = 4

# Compute Levenshtein edit distance between two sequences
def levenshtein(seq1, seq2):
    return Levenshtein.distance(seq1, seq2)

# Read DNA sequences from file
def read_sequences(filepath):
    with open(filepath, 'r') as file:
        sequences = file.read().splitlines()
    return sequences

# Manually compute pairwise distance matrix
def calculate_distance_matrix(sequences):
    num_sequences = len(sequences)
    distance_matrix = np.zeros((num_sequences, num_sequences))
    for i in range(num_sequences):
        for j in range(i + 1, num_sequences):
            dist = levenshtein(sequences[i], sequences[j])
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist
    return distance_matrix

# Option: Spectral Clustering (commented out)
# from sklearn.cluster import SpectralClustering
# from sklearn.metrics import pairwise_kernels
# def cluster_sequences(...): ...

# Use KMeans clustering
from sklearn.cluster import KMeans
def cluster_sequences(sequences, num_clusters):
    if len(sequences) < 2:
        return [0, 1] * len(sequences)
    else:
        distance_matrix = calculate_distance_matrix(sequences)
        feature_matrix = distance_matrix  # Raw distances used as features
        kmeans = KMeans(n_clusters=num_clusters, n_init=10, random_state=42)
        clustering = kmeans.fit(feature_matrix)
        return clustering.labels_

# Option: AgglomerativeClustering (commented out for alternative use)
# def cluster_sequences(...): ...

# Save first-level clustering result into individual files
def save_cluster_results_first(labels, sequences, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    clusters = {i: [] for i in range(max(labels) + 1)}
    for label, sequence in zip(labels, sequences):
        clusters[label].append(sequence)
    for cluster_id, seqs in clusters.items():
        with open(os.path.join(output_dir, f'cluster_{cluster_id}.txt'), 'w') as file:
            file.write('\n'.join(seqs))
    return clusters

# Save second-level binary clustering result into two folders
def save_cluster_results_second(labels, sequences, output_dir, cluster_id):
    sub_output_dir_0 = os.path.join(output_dir, 'cluster_double_result_0')
    sub_output_dir_1 = os.path.join(output_dir, 'cluster_double_result_1')
    os.makedirs(sub_output_dir_0, exist_ok=True)
    os.makedirs(sub_output_dir_1, exist_ok=True)

    clusters = {0: [], 1: []}
    if len(sequences) < 2:
        clusters[0].append(sequences[0])
        clusters[1].append(sequences[0])
    else:
        for label, sequence in zip(labels, sequences):
            clusters[label].append(sequence)

    if len(clusters[0]) < len(clusters[1]):
        clusters[0], clusters[1] = clusters[1], clusters[0]

    for sub_cluster_id in range(2):
        file_path = os.path.join(
            sub_output_dir_0 if sub_cluster_id == 0 else sub_output_dir_1,
            f'cluster_{cluster_id}_{sub_cluster_id}.txt'
        )
        with open(file_path, 'w') as file:
            file.writelines('\n'.join(clusters[sub_cluster_id]) + '\n')

    return clusters

# Save ground-truth sequence labels for purity evaluation
def separate_and_save_standard_sequences(filepath, output_dir):
    sequences = read_sequences(filepath)
    cluster1 = []
    cluster2 = []
    group_size = depth * n
    for i in range(0, len(sequences), group_size):
        cluster2.extend(sequences[i:i+depth])
        cluster1.extend(sequences[i+depth:i+group_size])
    with open(os.path.join(output_dir, 'standard_cluster1.txt'), 'w') as f:
        f.writelines('\n'.join(cluster1))
    with open(os.path.join(output_dir, 'standard_cluster2.txt'), 'w') as f:
        f.writelines('\n'.join(cluster2))
    print('Number of sequences in standard_cluster1:', len(cluster1))
    print('Number of sequences in standard_cluster2:', len(cluster2))
    return cluster1, cluster2

# Compute cluster purity
def calculate_purity(predicted_clusters, standard_cluster1, standard_cluster2):
    purity_stats = {}
    total_purity = 0
    for cluster_id, seqs in predicted_clusters.items():
        count_in_cluster1 = sum(seq in standard_cluster1 for seq in seqs)
        count_in_cluster2 = sum(seq in standard_cluster2 for seq in seqs)
        cluster_purity = max(count_in_cluster1, count_in_cluster2) / len(seqs) if len(seqs) > 0 else 0
        purity_stats[cluster_id] = {
            "total_sequences": len(seqs),
            "purity": cluster_purity,
            "details": {
                "count_in_standardcluster1": count_in_cluster1,
                "count_in_standardcluster2": count_in_cluster2
            }
        }
        total_purity += cluster_purity * len(seqs)
    total_sequences = sum(len(seqs) for seqs in predicted_clusters.values())
    overall_purity = total_purity / total_sequences if total_sequences > 0 else 0
    return overall_purity, purity_stats

# Step 1: Cluster all sequences into 10 groups
file_path = "...\\encoding\\errorstrands_output.txt"
output_directory = '...\\decoding\\clusterresult'
num_clusters = 10

sequences = read_sequences(file_path)
labels = cluster_sequences(sequences, num_clusters)
clusters = save_cluster_results_first(labels, sequences, output_directory)

for cluster_id, seqs in clusters.items():
    print(f"Cluster {cluster_id}: {len(seqs)} sequences")

print("Initial clustering saved (divided into 10 groups).")

# Optional: visualize per-cluster purity
def plot_purity_barchart(purity_stats):
    cluster_ids = list(purity_stats.keys())
    purity_values = [purity_stats[cluster_id]["purity"] for cluster_id in cluster_ids]
    plt.bar(cluster_ids, purity_values, color='skyblue')
    plt.xlabel('Cluster ID')
    plt.ylabel('Purity')
    plt.title('Purity of Each Cluster')
    plt.ylim(0, 1)
    plt.xticks(cluster_ids)
    plt.show()

# Step 2: Perform binary clustering on each group and evaluate purity
def main():
    error_strands_file = '...\\encoding\\errorstrands_output.txt'
    standard_output_dir = '...\\decoding'
    input_directory = '...\\decoding\\clusterresult'
    num_clusters = 2

    standard_cluster1, standard_cluster2 = separate_and_save_standard_sequences(error_strands_file, standard_output_dir)

    total_purity_list = []
    purity_stats = {}

    for cluster_id in range(10):
        input_file = os.path.join(input_directory, f'cluster_{cluster_id}.txt')
        sequences = read_sequences(input_file)
        labels = cluster_sequences(sequences, num_clusters)
        clusters = save_cluster_results_second(labels, sequences, standard_output_dir, cluster_id)

        overall_purity, stats = calculate_purity(clusters, standard_cluster1, standard_cluster2)
        total_purity_list.append(overall_purity)
        purity_stats[cluster_id] = {"purity": overall_purity}

        print(f"Cluster {cluster_id} overall purity: {overall_purity:.2%}")
        for k, v in stats.items():
            print(f"Cluster {cluster_id}_{k}: {v['total_sequences']} sequences, purity: {v['purity']:.2%}, "
                  f"details: {v['details']}")

    total_mean_purity = np.mean(total_purity_list)
    print(f"Total purity: {total_mean_purity:.2%}")

    # plot_purity_barchart(purity_stats)

if __name__ == "__main__":
    main()
