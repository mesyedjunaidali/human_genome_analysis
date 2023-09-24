import pandas as pd
import matplotlib.pyplot as plt

protein_sequences = {
    'Alanine': ['GCU', 'GCC', 'GCA', 'GCG'],
    'Arginine': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'Asparagine': ['AAU', 'AAC'],
    'Aspartic Acid': ['GAU', 'GAC'],
    'Cysteine': ['UGU', 'UGC'],
    'Glutamic Acid': ['GAA', 'GAG'],
    'Glutamine': ['CAA', 'CAG'],
    'Glycine': ['GGU', 'GGC', 'GGA', 'GGG'],
    'Histidine': ['CAU', 'CAC'],
    'Isoleucine': ['AUU', 'AUC', 'AUA'],
    'Leucine': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'Lysine': ['AAA', 'AAG'],
    'Methionine': ['AUG'],
    'Phenylalanine': ['UUU', 'UUC'],
    'Proline': ['CCU', 'CCC', 'CCA', 'CCG'],
    'Serine': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    'Threonine': ['ACU', 'ACC', 'ACA', 'ACG'],
    'Valine': ['GUU', 'GUC', 'GUA', 'GUG']
}

# Initialize counters for each protein sequence
protein_counts = {protein: 0 for protein in protein_sequences.keys()}

# Initialize counters for A+T and G+C
a_t_count = 0
g_c_count = 0

# Initialize counters for individual nucleotides
nucleotide_counts = {'A': 0, 'G': 0, 'T': 0, 'C': 0}

# Read the large text file line by line
with open('GRCh38_latest_genomic_sequence.txt', 'r') as file:
    for line in file:
        sequence = line.strip()
        sequence = sequence.upper()  # Convert to uppercase for case-insensitive matching

        # Update protein counts
        for protein, sequences in protein_sequences.items():
            for seq in sequences:
                protein_counts[protein] += sequence.count(seq)

        # Update A+T and G+C counts
        a_t_count += sequence.count('A') + sequence.count('T')
        g_c_count += sequence.count('G') + sequence.count('C')

        # Update nucleotide counts
        for nucleotide in sequence:
            if nucleotide in nucleotide_counts:
                nucleotide_counts[nucleotide] += 1

# Calculate the total count and percentage of each protein
total_count = sum(protein_counts.values())
percentage_compositions = {protein: (count / total_count * 100) for protein, count in protein_counts.items()}

# Calculate the (A+T)/(G+C) ratio
at_gc_ratio = a_t_count / g_c_count

# Calculate the total count of nucleotides
total_nucleotides = sum(nucleotide_counts.values())

# Calculate the percentage composition of individual nucleotides
percentage_nucleotides = {nucleotide: (count / total_nucleotides * 100) for nucleotide, count in nucleotide_counts.items()}

# Convert nucleotide compositions to a pandas DataFrame
nucleotide_df = pd.DataFrame.from_dict(percentage_nucleotides, orient='index', columns=['Percentage'])
nucleotide_df = nucleotide_df.sort_index()

# Convert protein compositions to a pandas DataFrame
protein_df = pd.DataFrame.from_dict(percentage_compositions, orient='index', columns=['Percentage'])
protein_df = protein_df.sort_values(by='Percentage', ascending=False)

# Plotting nucleotide composition
nucleotide_df.plot(kind='line', marker='o')
plt.xlabel('Nucleotide')
plt.ylabel('Percentage Composition')
plt.title('Nucleotide Composition')
plt.show()

# Plotting protein composition
protein_df.plot(kind='line', marker='o')
plt.xlabel('Protein')
plt.ylabel('Percentage Composition')
plt.title('Protein Composition')
plt.show()
