import os
import numpy as np
import Levenshtein
from sklearn.cluster import AgglomerativeClustering

depth = 60
n = 5

# Compute edit distance using Levenshtein library
def levenshtein(seq1, seq2):
    return Levenshtein.distance(seq1, seq2)

# Read DNA sequences from a text file
def read_sequences(filepath):
    with open(filepath, 'r') as file:
        sequences = file.read().splitlines()
    return sequences

# Manually compute the pairwise distance matrix
def calculate_distance_matrix(sequences):
    num_sequences = len(sequences)
    distance_matrix = np.zeros((num_sequences, num_sequences))
    for i in range(num_sequences):
        for j in range(i + 1, num_sequences):
            dist = levenshtein(sequences[i], sequences[j])
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist
    return distance_matrix

# Perform agglomerative clustering
def cluster_sequences(sequences, num_clusters):
    distance_matrix = calculate_distance_matrix(sequences)
    clustering = AgglomerativeClustering(n_clusters=num_clusters, metric='precomputed', linkage='average')
    clustering.fit(distance_matrix)
    return clustering.labels_

# Save clustering results to separate files
def save_cluster_results(labels, sequences, output_dir, cluster_id):
    clusters = {i: [] for i in range(max(labels) + 1)}
    for label, sequence in zip(labels, sequences):
        clusters[label].append(sequence)
    for sub_cluster_id, seqs in clusters.items():
        sub_cluster_filename = f'cluster_{cluster_id}_{sub_cluster_id}.txt'
        with open(os.path.join(output_dir, sub_cluster_filename), 'w') as file:
            file.write('\n'.join(seqs))
    return clusters

# Split and save standard sequences for comparison
def separate_and_save_standard_sequences(filepath, output_dir):
    sequences = read_sequences(filepath)
    cluster1 = []
    cluster2 = []
    group_size = depth * n  # Each group contains 60 normal + 4*60 misleading sequences (if n=4)
    for i in range(0, len(sequences), group_size):
        cluster2.extend(sequences[i:i+depth])  # First 60 sequences are normal (assigned to cluster2)
        cluster1.extend(sequences[i+depth:i+group_size])  # Remaining sequences are misleading (cluster1)
    with open(os.path.join(output_dir, 'standard_cluster1.txt'), 'w') as f:
        f.writelines('\n'.join(cluster1))
    with open(os.path.join(output_dir, 'standard_cluster2.txt'), 'w') as f:
        f.writelines('\n'.join(cluster2))
    print('Number of sequences in standard_cluster1:', len(cluster1))
    print('Number of sequences in standard_cluster2:', len(cluster2))
    return cluster1, cluster2

# Calculate clustering purity
def calculate_purity(predicted_clusters, standard_cluster1, standard_cluster2):
    purity_stats = {}
    total_purity = 0
    for cluster_id, seqs in predicted_clusters.items():
        count_in_cluster1 = sum(seq in standard_cluster1 for seq in seqs)
        count_in_cluster2 = sum(seq in standard_cluster2 for seq in seqs)
        # Purity is defined as the proportion of majority-class sequences in each cluster
        cluster_purity = max(count_in_cluster1, count_in_cluster2) / len(seqs) if len(seqs) > 0 else 0
        purity_stats[cluster_id] = {
            "total_sequences": len(seqs),
            "purity": cluster_purity,
            "details": {
                "count_in_standardcluster1": count_in_cluster1,
                "count_in_standardcluster2": count_in_cluster2
            }
        }
        total_purity += cluster_purity * len(seqs)  # Weighted purity

    total_sequences = sum(len(seqs) for seqs in predicted_clusters.values())
    overall_purity = total_purity / total_sequences if total_sequences > 0 else 0
    return overall_purity, purity_stats

# Main program
error_strands_file = '...\\encoding\\errorstrands_output.txt'
standard_output_dir = '...\\decoding'
input_directory = '...\\decoding\\clusterresult'
output_directory = '...\\decoding\\clusteranalysis'
num_clusters = 2  # Cluster sequences in each file into 2 groups

# Split and save standard clusters
standard_cluster1, standard_cluster2 = separate_and_save_standard_sequences(error_strands_file, standard_output_dir)

# Analyze cluster files and compute purity
total_purity_list = []
for cluster_id in range(10):
    input_file = os.path.join(input_directory, f'cluster_{cluster_id}.txt')
    sequences = read_sequences(input_file)
    labels = cluster_sequences(sequences, num_clusters)
    clusters = save_cluster_results(labels, sequences, output_directory, cluster_id)

    # Compute purity and output detailed statistics
    overall_purity, stats = calculate_purity(clusters, standard_cluster1, standard_cluster2)
    total_purity_list.append(overall_purity)
    print(f"Cluster {cluster_id} overall purity: {overall_purity:.2%}")
    for k, v in stats.items():
        print(f"Cluster {cluster_id}_{k}: {v['total_sequences']} sequences, purity: {v['purity']:.2%}, "
              f"details: {v['details']}")

# Compute average purity
total_mean_purity = np.mean(total_purity_list)
print(f"Total purity: {total_mean_purity:.2%}")
