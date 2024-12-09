#with open("./SRR000002.fasta", "r") as fasta, open("./input.tsv", "w") as tsv:
#    for line in fasta:
#        if not line.startswith(">"):
#            tsv.write(line.strip() + "\n")

with open("./SRR000002.fasta", "r") as infile, open("./input.tsv", "w") as outfile:
    sequence = []
    for line in infile:
        if line.startswith(">"):  # Ligne d'en-tête
            if sequence:
                outfile.write("".join(sequence) + "\n")
                sequence = []  # Réinitialiser la séquence
        else:
            sequence.append(line.strip())  # Ajouter la séquence sans saut de ligne
    if sequence:  # Ajouter la dernière séquence
        outfile.write("".join(sequence) + "\n")
