# 2024_Meta16S_minion_Morgane

## Summary
[TOC]

## Emplacement données brutes
Les données brutes sont stockées sur le QNAP2. Les flowcells sont :
1ere analyse:
- flowcell16s_240605
- flowcell16s_240606_2

2d analyse :
- flowcell16s_241022
- flowcell16s_241023
- flowcell16s_241028
- flowcell16s_241029



## 1ere analyse 08/07/2024 : flowcell16s_240605 & flowcell16s_240606_2
clone epi2me github repository (a nextflow piepline for long read 16S metagenomic using kraken2 or minimap2)
https://github.com/epi2me-labs/wf-metagenomics

### database
Utilisation de la database k2_standard_08gb_20231009.tar.gz (5,5Go)


Lancement sur les 48 échantillons T0 à T2:



### Lancement Epi2me

Lancé avec la base de données k2_standard_08gb_20231009.tar.gz (database) et new_taxdump_11062024.tar.gz (taxonomy).

Avant de lancer le pipeline, il faut que chaque échantillons soit placé dans le dossier Fastq et dans un dossier à son nom pour bien etre pris en compte par le pipeline.

```bash
cd /mnt/SSD_02/12_2024_Meta16S_minion_Morgane/
nextflow run wf-metagenomics/ --fastq Fastq/ -profile singularity
```


### R analysis

Il faut supprimer les "'" qui pose probleme dans le script R
```bash
sed "s/'//g" abundance_table_species.tsv > abundance_table_species_clean.tsv
```


=> Voir script **phyloseq.R**

#### Résultats

Les résutlats / figures sont enregistrées dans le ppt : Présentation_Projet.pptx
K:\Personnel_non_permanent\Thibault\12_metagenomique16_minion_Morgane





## 2eme analyse 11/07/2024 => Test database 60Go : FAIL

https://benlangmead.github.io/aws-indexes/k2

A partir du lien ci-dessus, téléchargement de la base de donnée 60Go => voir si cela change quelque chose à nos resultats
NE FONCTIONNE PAS => Plante systématiquement




## 3eme analyse 14/11/2024 : 6 flowcells

### Formatage des donnnées

#### concatenation des fastq.gz en 1 seul fichier + rename
Les nouvelles données sont réparties dans 4 dossiers correspondant au 4 flowcells. Pour chacune des flowcell, nous avons 1 dossier par barcodes et dans ce dossier plusieurs fichiers ".fastq.gz" qu'il faut concatener. Pour optimiser la concaténation de ces fichiers, nous allons créer un script permettant d'avoir 1 fichier final par barcode qui a pour nom "Barcode_NumRun.fastq.gz"

Copier ce script dans le dossier désiré, remplacer le PATH de "base_dir", changer les permissions avec chmod +x et le lancer.

```bash
#!/bin/bash

# Remplacer par le chemin du répertoire contenant les dossiers de barcode
base_dir="PATH"

# Parcourir chaque dossier commençant par "barcode"
for barcode_dir in "$base_dir"/barcode*/; do
    # Récupérer le nom du dossier
    barcode=$(basename "$barcode_dir")

    # Extraire le numéro du run depuis le premier fichier du dossier
    run_id=$(ls "$barcode_dir"/*.fastq.gz | head -n 1 | sed -n 's/.*_runid_\([a-z0-9]*\)_.*/\1/p')

    # Créer le nom de fichier final
    output_file="${base_dir}/${barcode}_${run_id}.fastq.gz"

    # Concaténer tous les fichiers fastq.gz dans le fichier final
    cat "$barcode_dir"/*.fastq.gz > "$output_file"

    echo "Fichier concaténé pour ${barcode_dir}: ${output_file}"
done
```


Il faut ensuite renommer ces fichiers avec les metadata (faire un fichier de correspondance et renommer automatiquement avec un "mv" ?)


#### Import des données sur wlyo3

```bash
# Pour les 4 premieres flowcells
cd /mnt/NGS/longread
nohup rsync flowcell16s_241022/basecalling/pass/concat_fastq/* /mnt/SSD_02/12_2024_Meta16S_minion_Morgane/Fastq/flowcell22/ &
nohup rsync flowcell16s_241023/basecalling/pass/concat_fastq/* /mnt/SSD_02/12_2024_Meta16S_minion_Morgane/Fastq/flowcell23/ &
nohup rsync flowcell16s_241028/basecalling/pass/concat_fastq/* /mnt/SSD_02/12_2024_Meta16S_minion_Morgane/Fastq/flowcell28/ &
nohup rsync flowcell16s_241029/basecalling/pass/concat_fastq/* /mnt/SSD_02/12_2024_Meta16S_minion_Morgane/Fastq/flowcell29/ &

# Pour les 2 premieres flowcells
nohup rsync flowcell16s_240605_r10.4.1_SQK-16S114-24_dorado0.4.2_e8.2_400bps_sup_500/dorado_barcoder/*fastq.gz /mnt/SSD_02/12_2024_Meta16S_minion_Morgane/Fastq/flowcell05/ &
nohup rsync flowcell16s_240606_2_r10.4.1_SQK-16S114-24_dorado0.4.2_e8.2_400bps_sup_500/dorado_barcoder/*fastq.gz /mnt/SSD_02/12_2024_Meta16S_minion_Morgane/Fastq/flowcell06_2/ &
```

#### Formater les Fastq pour Epi2me

Certains échantillons sont non appariés : on les retire de l'analyse. D'autres on été séquencé 2 fois : on ajoute un "_bis" au noms de l'échantillon.

Morgane a envoyé un fichier avec les échantillon à supprimer/modifier:

Nom_tubes       Flowcell        Action
T1_11   1       à supprimer
T1_12   1       à supprimer
T1_13   1       à supprimer
T1_14   1       à supprimer
T1_15   1       à supprimer
T1_16   2       à supprimer
T1_17   2       à supprimer
T1_18   2       à supprimer
T1_19   2       à supprimer
T1_20   2       à supprimer
T1_21   2       à supprimer
T2_11   2       à supprimer
T2_12   2       à supprimer
T2_13   2       à supprimer
T2_14   2       à supprimer
T2_15   2       à supprimer
T2_16   2       à supprimer
T2_17   2       à supprimer
T2_18   2       à supprimer
T2_5    1
T2_7    1
T2_5    3       "Ajouter bis
échantillon avec moins de 3000 reads"
T2_7    3       Ajouter bis
T2_22   5
T2_23   5
T2_22   6       Ajouter bis
T2_23   6       Ajouter bis


```bash
# Ajout un "bis" au nom de l'échantillon juste apres le chiffre pour pouvoir récuperer le pattern dans la commande ci-dessous
# De la meme manière rennomemr les Controle positif
# il y en a que 4 dont on doit modifier le nom.

# Pour chaque fichier fastq.gz dans les dossier flowcell, on modifie le nom pour ne garder que les 2 premier pattern avant le second "_". On les copie ensuite dans un dossier portant le nom de l'échantillon.
for flowcell_dir in flowcell*; do
    cd "$flowcell_dir"
    for file in *.fastq.gz; do
        sample_name=$(echo "$file" | cut -d'_' -f1,2)
        mkdir -p "../$sample_name"
        cp "$file" "../$sample_name/"
    done
    cd ..
done

# Il devrait y avoir 124 samples
# il faut déplacer les dossier flowcell dans un autre dossier sinon cela plante lors du lancement du pipeline
```


### Lancement Epi2me

Lancé avec la base de données k2_standard_08gb_20231009.tar.gz (database) et new_taxdump_11062024.tar.gz (taxonomy).

Avant de lancer le pipeline, il faut que chaque échantillons soit placé dans le dossier Fastq et dans un dossier à son nom pour bien etre pris en compte par le pipeline => voir étapes réalisées ci-dessus.

```bash
cd /mnt/SSD_02/12_2024_Meta16S_minion_Morgane/
nextflow run wf-metagenomics/ --fastq Fastq/ -profile singularity

cd output
sed "s/'//g" abundance_table_species.tsv > abundance_table_species_clean.csv
```

### R analysis

#### Résultats
