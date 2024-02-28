from collections import Counter
import numpy as np

def read_fasta(arq):
    seq = ''
    with open(arq) as f:
        f.readline()  # Ignora a primeira linha (cabeçalho)
        for line in f:
            seq += line.strip()  # Adiciona cada linha da sequência, removendo espaços em branco
    return seq

def comparacao_simples(seq1, seq2):
    matches = sum(a == b for a, b in zip(seq1, seq2))
    return matches

def contar_aminoacidos(seq):
    return Counter(seq)

def distancia_manhattan(v1, v2):
    return np.sum(np.abs(v1 - v2))

def distancia_euclidiana(v1, v2):
    return np.sqrt(np.sum((v1 - v2) ** 2))

def distancia_supremum(v1, v2):
    return np.max(np.abs(v1 - v2))

def similaridade_cosseno(v1, v2):
    return np.dot(v1, v2) / (np.sqrt(np.dot(v1, v1)) * np.sqrt(np.dot(v2, v2)))

def preparar_vetores(contagem1, contagem2):
    todos_aminoacidos = set(contagem1.keys()) | set(contagem2.keys())
    v1 = np.array([contagem1.get(amino, 0) for amino in todos_aminoacidos])
    v2 = np.array([contagem2.get(amino, 0) for amino in todos_aminoacidos])
    return v1, v2

# Leitura das sequências dos arquivos FASTA
seq_hamster = read_fasta('hamster.fasta')
seq_rat = read_fasta('rat.fasta')
seq_horse = read_fasta('horse.fasta')

# Pares de organismos para comparação
pares = [('Hamster', seq_hamster, 'Rato', seq_rat),
         ('Hamster', seq_hamster, 'Cavalo', seq_horse),
         ('Rato', seq_rat, 'Cavalo', seq_horse)]

# Executar comparações
for org1_nome, seq1, org2_nome, seq2 in pares:
    print(f"\nComparando {org1_nome} VS {org2_nome}:")
    
    # Comparação simples
    matches = comparacao_simples(seq1, seq2)
    print(f"  Matches Simples: {matches}")
    
    # Contagem de aminoácidos/nucleotídeos
    contagem1 = contar_aminoacidos(seq1)
    contagem2 = contar_aminoacidos(seq2)
    
    # Preparar vetores
    v1, v2 = preparar_vetores(contagem1, contagem2)
    
    # Calcular distâncias e similaridade
    manhattan = distancia_manhattan(v1, v2)
    euclidiana = distancia_euclidiana(v1, v2)
    supremum = distancia_supremum(v1, v2)
    cosseno = similaridade_cosseno(v1, v2)
    
    # Imprimir resultados
    print(f"  Distância Manhattan: {manhattan}")
    print(f"  Distância Euclidiana: {euclidiana}")
    print(f"  Distância Supremum: {supremum}")
    print(f"  Similaridade de Cosseno: {cosseno}")
