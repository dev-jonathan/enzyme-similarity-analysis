from collections import Counter
import numpy as np

# Função para ler arquivos FASTA
def read_fasta(arq):
    seq = ''
    with open(arq) as f:
        f.readline()  # Ignora a primeira linha (cabeçalho)
        for line in f:
            seq += line.strip()  # Adiciona cada linha da sequência, removendo espaços em branco
    return seq

# Função para comparação simples
def comparacao_simples(seq1, seq2):
    matches = sum(a == b for a, b in zip(seq1, seq2))
    return matches

# Função para contar ocorrências de aminoácidos (ou nucleotídeos, neste caso)
def contar_aminoacidos(seq):
    return Counter(seq)

# Funções para cálculo de distâncias e similaridade
def distancia_manhattan(v1, v2):
    return np.sum(np.abs(v1 - v2))

def distancia_euclidiana(v1, v2):
    return np.sqrt(np.sum((v1 - v2) ** 2))

def distancia_supremum(v1, v2):
    return np.max(np.abs(v1 - v2))

def similaridade_cosseno(v1, v2):
    return np.dot(v1, v2) / (np.sqrt(np.dot(v1, v1)) * np.sqrt(np.dot(v2, v2)))

# Função auxiliar para preparar vetores numéricos
def preparar_vetores(contagem1, contagem2):
    todos_nucleotideos = set(contagem1.keys()) | set(contagem2.keys())
    v1 = np.array([contagem1.get(nucleotideo, 0) for nucleotideo in todos_nucleotideos])
    v2 = np.array([contagem2.get(nucleotideo, 0) for nucleotideo in todos_nucleotideos])
    return v1, v2

# Leitura das sequências dos arquivos FASTA
seq_hamster = read_fasta('hamster.fasta')
seq_rat = read_fasta('rat.fasta')
seq_horse = read_fasta('horse.fasta')

# Lista de organismos e suas sequências para facilitar a iteração
organismos = [('Hamster', seq_hamster), ('Rato', seq_rat), ('Cavalo', seq_horse)]

# Comparação de proximidade e cálculo de distâncias/similaridades
for i in range(len(organismos)):
    for j in range(i + 1, len(organismos)):
        org1_name, seq1 = organismos[i]
        org2_name, seq2 = organismos[j]
        
        # Comparação simples
        simples_matches = comparacao_simples(seq1, seq2)
        
        # Contagem de ocorrências
        contagem1 = contar_aminoacidos(seq1)
        contagem2 = contar_aminoacidos(seq2)
        
        # Converter contagens em vetores para cálculo de distâncias/similaridades
        v1, v2 = preparar_vetores(contagem1, contagem2)
        
        # Calcula as distâncias e similaridade
        manhattan = distancia_manhattan(v1, v2)
        euclidiana = distancia_euclidiana(v1, v2)
        supremum = distancia_supremum(v1, v2)
        cosseno = similaridade_cosseno(v1, v2)
        
        # Imprime os resultados
        print(f"\nComparação {org1_name} vs {org2_name}:")
        print(f"Matches Simples: {simples_matches}")
        print(f"Distância Manhattan: {manhattan}")
        print(f"Distância Euclidiana: {euclidiana}")
        print(f"Distância Supremum: {supremum}")
        print(f"Similaridade de Cosseno: {cosseno}")
