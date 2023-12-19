from Bio import SeqIO
import numpy as np
from newick import loads

def read_sequences(dataset_number):
    """
    Read sequences from FASTA files in the specified dataset folder.

    Parameters:
    - dataset_number (int): The number of the dataset.

    Returns:
    - sequences (dict): A dictionary containing sequence names as keys and sequences as values.
    """
    folder_name = f"dataset_{dataset_number}"

    sequences = {}
    for seq_number in range(1, 11):
        file_name = f"seq_{seq_number}.fasta"
        file_path = f"{folder_name}/{file_name}"

        try:
            record = SeqIO.read(file_path, "fasta")
            sequences[file_name] = str(record.seq)[1:]
        except FileNotFoundError:
            print(f"Arquivo {file_name} não encontrado.")

    return sequences

def needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-2):
    """
    Perform Needleman-Wunsch algorithm for sequence alignment.

    Parameters:
    - seq1 (str): First sequence.
    - seq2 (str): Second sequence.
    - match_score (int): Score for a match.
    - mismatch_score (int): Score for a mismatch.
    - gap_penalty (int): Penalty for a gap.

    Returns:
    - int: Alignment score.
    """
    rows, cols = len(seq1) + 1, len(seq2) + 1
    matrix = np.zeros((rows, cols), dtype=int)

    for i in range(1, rows):
        matrix[i][0] = i * gap_penalty
    for j in range(1, cols):
        matrix[0][j] = j * gap_penalty

    for i in range(1, rows):
        for j in range(1, cols):
            match = matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score)
            delete = matrix[i - 1][j] + gap_penalty
            insert = matrix[i][j - 1] + gap_penalty
            matrix[i][j] = max(match, delete, insert)

    return matrix[-1][-1]

def calculate_alignment_scores(sequences, limit=10000):
    """
    Calculate alignment scores for all pairs of sequences.

    Parameters:
    - sequences (dict): Dictionary of sequences.
    - limit (int): Limit for sequence length.

    Returns:
    - list: List of alignment scores.
    """
    scores = []

    for i in range(1, 11):
        for j in range(i + 1, 11):
            seq1 = sequences[f"seq_{i}.fasta"][:limit]
            seq2 = sequences[f"seq_{j}.fasta"][:limit]
            score = needleman_wunsch(seq1, seq2)
            scores.append(f"seq_{i} vs seq_{j}: {score}")

    return scores

def read_distance_matrix(dataset_number):
    """
    Read distance matrix from a file.

    Parameters:
    - dataset_number (int): The number of the dataset.

    Returns:
    - otu_labels (list): List of labels for OTUs.
    - distance_matrix (numpy.ndarray): Distance matrix.
    """
    file_path = f"dataset_{dataset_number}/dist_matriz.txt"
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            otu_labels = lines[0].split()[1:]
            otu_labels.insert(0, 'seq_1')
            distance_matrix = []

            for line in lines[1:]:
                tokens = line.split()
                row_values = list(map(int, tokens[1:]))
                distance_matrix.append(row_values)

            return otu_labels, np.array(distance_matrix)

    except FileNotFoundError:
        print(f"O arquivo {file_path} não foi encontrado.")
        return None

def upgma(distance_matrix, labels):
    """
    Perform UPGMA clustering on a distance matrix.

    Parameters:
    - distance_matrix (numpy.ndarray): Distance matrix.
    - labels (list): List of labels for the data points.

    Returns:
    - str: Newick format tree string.
    """
    while distance_matrix.shape[0] > 1:
        n = distance_matrix.shape[0]
        i, j = find_lowest_value(distance_matrix)

        dist_node_to_u = distance_matrix[i, j] / 2
        new_node = f"({labels[i]}: {round(dist_node_to_u, 4)}, {labels[j]}: {round(dist_node_to_u, 4)})"
        new_labels = [new_node]

        for label in labels:
            if label not in {labels[i], labels[j]}:
                new_labels.append(label)

        dprime = np.zeros((n - 1, n - 1))
        ij_indexes = [i, j]
        dprime[1:, 1:] = np.delete(np.delete(distance_matrix, ij_indexes, axis=1),
                                   ij_indexes, axis=0)

        dist_k_to_u = (distance_matrix[i] + distance_matrix[j]) / 2
        dist_k_to_u = np.delete(dist_k_to_u, ij_indexes)
        dist_k_to_u = np.concatenate([[0], dist_k_to_u])

        for k in range(n - 1):
            dprime[0, k] = dprime[k, 0] = dist_k_to_u[k]

        distance_matrix = np.copy(dprime)
        labels = new_labels.copy()

    tree = labels[0]
    return tree + ";"

def find_lowest_value(matrix):
    """
    Find the indices of the lowest value in a 2D matrix.

    Parameters:
    - matrix (numpy.ndarray): 2D matrix.

    Returns:
    - tuple: Indices of the lowest value.
    """
    min_val_idx = np.unravel_index(np.argmin(matrix), matrix.shape)
    return min_val_idx

# Exemplo de uso
dataset_number = input("Enter the dataset number: ")
question = input("Enter the question (b or d): ")
# *question b
#   needs dataset_{number} folder
if question == 'b':
    sequences = read_sequences(int(dataset_number))
    alignment_scores = calculate_alignment_scores(sequences)
    for score in alignment_scores:
        print(score)
# *question d
#   needs in dataset_{number} dist_matriz.txt
if question == 'd':
    otu_labels, distance_matrix = read_distance_matrix(int(dataset_number))
    tree = upgma(distance_matrix, otu_labels)
    print(tree)
    print(loads(tree)[0].ascii_art())
