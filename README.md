# Assembleur basé sur les graphes de Debruijn

Vous trouverez la description complète du TP [ici]( 
https://docs.google.com/document/d/1P4v3bHbSurD7RXA-ldVwmtNKGvWsnBae51RMGye_KLs/edit?usp=sharing).

## Introduction

L’objectif de ce TP est d’assembler le génome de l’entérovirus A71. Ce génome présente l’intérêt d’être très court: 7408 nucléotides, linéaire et non segmenté.
Le fichier fastq dont vous disposez a été généré à l’aide du programme ART [Huang 2011] via la commande:
```bash
art_illumina -i eva71.fna -ef -l 100 -f 20 -o eva71 -ir 0 -dr 0 -ir2 0 -dr2 0 -na -qL 41 -rs 1539952693
```
Les lectures ont une qualité maximale (41) et ne présentent pas d’insertion. Seuls les lectures correspondant aux brins 5’ -> 3’ vous sont ici fournies.

## Fonctionnement du programme

L’assembleur suit les étapes classiques d’un **assemblage de novo** basé sur les graphes de De Bruijn :

1. Lecture des séquences depuis le fichier FASTQ.  
2. Découpage des lectures en *k*-mers.  
3. Construction du graphe de De Bruijn (`networkx.DiGraph`).  
4. Simplification du graphe (reconstruction des **unitigs**).  
5. Extraction et écriture des **contigs** finaux au format FASTA.  

## Données

Dans le dossier [`data/`](./data/), vous trouverez :

| Fichier | Description |
|----------|-------------|
| `eva71.fna` | Génome de référence de l’entérovirus A71 |
| `eva71_plus_perfect.fq` | Lectures simulées au format FASTQ |


## Installation des dépendances

Ce projet utilise [**Pixi**](https://pixi.sh) comme gestionnaire d’environnement et de dépendances. Toutes les dépendances sont listées dans [`pixi.toml`](./pixi.toml).

### Installation et mise à jour de l’environnement :

```bash
pixi install   # Installe les dépendances si nécessaire
pixi update    # Met à jour selon pixi.toml
```
Les principales librairies utilisées sont : :
- `networkx` — manipulation du graphe de De Bruijn  
- `pytest`, `pytest-cov` — tests unitaires et couverture  
- `pylint` — vérification syntaxique

## Utilisation

Le script principal est [`debruijn.py`](./debruijn/debruijn.py). Il prend en entrée un fichier FASTQ et assemble les **contigs** à partir des **k-mers**.

### Commande de base :

```bash
pixi run python debruijn/debruijn.py -i data/eva71_plus_perfect.fq -k 21 -o results/contigs.fasta
```
### Arguments disponibles

- `-h, --help`         : afficher l’aide et quitter
- `-i FASTQ_FILE`      : fichier FASTQ single-end (obligatoire)
- `-k KMER_SIZE`       : taille des k-mers (par défaut : 21)
- `-o OUTPUT_FILE`     : fichier de sortie contenant les contigs au format FASTA

---

## Structure du projet

```
├── README.md
├── data/
│   ├── eva71.fna
│   └── eva71_plus_perfect.fq
├── debruijn/
│   └── debruijn.py
├── pixi.lock
├── pixi.toml
└── tests/
```
## Résultats attendus

L’exécution du programme sur les données du TP génère un fichier `contigs.fasta` contenant les séquences assemblées.  
Ces contigs peuvent ensuite être comparés :

- au **génome de référence** [`data/eva71.fna`](./data/eva71.fna),  
- ou à d’autres outils d’assemblage comme **SPAdes** ou **Velvet** pour évaluer la qualité de la reconstruction.

## Outils complémentaires

Pour les analyses de validation, la suite **BLAST+** peut être utilisée afin de comparer les contigs au génome de référence :  

```bash
makeblastdb -in data/eva71.fna -dbtype nucl -out data/eva71_db
blastn -query results/contigs.fasta -db data/eva71_db -out results/blast_results.txt
```

## Contact

En cas de questions, vous pouvez me contacter par email: [assa.diabira@etu.u-paris.fr](mailto:assa.diabira@etu.u-paris.fr)
