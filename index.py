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


# Abrindo o arquivo log_resul.txt para escrita
with open('log_resultado.txt', 'w') as log_file:
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

            # Imprime os resultados no arquivo log_resul.txt
            print(f"===========================================", file=log_file)
            print(f"|Comparação {org1_name} vs {org2_name}:|", file=log_file)
            print(f"Matches Simples: {simples_matches}", file=log_file)
            print(f"Distância Manhattan: {manhattan}", file=log_file)
            print(f"Distância Euclidiana: {euclidiana}", file=log_file)
            print(f"Distância Supremum: {supremum}", file=log_file)
            print(f"Similaridade de Cosseno: {cosseno}", file=log_file)