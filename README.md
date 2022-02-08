# The Great Cattle Antibody Repertoire
Longtitudal study of vaccination of 204 calves against the bovine respiratory disease 

## Repository files
- **aux:** antibody titers, calf ages, V gene usages, ultralong CDRH3 fractions;
- **gsv_finding:** intermediate results of IGHV gene genotyping, detected germline and somatic variations (GSVs) and their R-ratios;
- **igqtl_results:** clusters of calves computed using R-ratios of the detected GSVs.

## AIRR metadata
Metadata for 204 processed antibody repertoires are available [here](https://livejohnshopkins-my.sharepoint.com/:f:/g/personal/isafono1_jh_edu/Egd-hagOhZ1Nq_0KRrmlYgkBdf5WJ048SXeI5JvrJQFgGQ?e=dWefCV). 

Each zipped CSV annotation file describes four Rep-Seq libraries corresponding to a single animal and includes all distinct VDJ nucleotide sequences. Each line represents a distinct nucleotide sequence and contains with the following information:

- Nucleotide and amino acid sequences,
- Read count,
- Time point (–3, 0, 3, or 6),
- IGHV and IGHJ gene hits,
- Nucleotide and amino acid sequences of CDR1, CDR2, CDR3,
- Whether or not CDR3 is ultralong (i.e., the length of CDR3 exceeds 150 nt),
- Clone ID. Clones are specific to time point, i.e., they are computed for sequences from a given time point. A clone is defined as a group of sequences sharing V hit and J hit and having identical CDR3s.
- Alignment of the variable gene segment against the closest germline IGHV gene.

## Germline and somatic variation data
- **vgene_genotyping:** the directory contains subdirectories corresponding to 204 cattle subjects. Each subdirectory consists of TXT files with nucleotide counts for all positions of all V genes. Each file corresponds to a single V gene and contains a table showing the numbers of clones (in this case, clones = sequences sharing V and J genes and having identical CDR3s) with A, C, G, and T nucleotides at positions 0,...,L-1, where L is the length of the V gene. A fragment of the 14007_mismatches/IGHV1-7_mismatches.txt file is shown below. E.g., 4, 7, 10010, and 10 clones have A, C, G, and T at position 38 (0-based) of IGHV1-7, respectively.
 
| Position | NumAs | NumCs | NumGs | NumTs |
| --- | --- | --- | --- | --- |
| ... |||||
|38 | 4 | 7 | 10010 | 10 |
|39 | 1 | 9978 | 2 | 50 |
|40 | 0 | 10018 | 6 | 7 |
|41 | 24 | 3664 | 6332 | 11 |
|42 | 5 | 24 | 7 | 9995 |
|...|||||

- **gsv_r_ratios.txt**: TXT files containing R-ratios for 52 GSVs detected by IgQTL tool using the nucleotide counts from **vgene_genotyping**. 

## Scripts
- **igqtl.py**: takes **vgene_genotyping** directory as an input, computes germline and somatic variations, and reports their R-ratios for all subjects:
```
python igqtl.py vgene_genotyping_dir output_dir
```

Output directory contains R-ratios of GSVs per V gene and all GSV combined together. Each identified GSV is written as a column and have following format:
```
VGENE:POSITION_N1N2 
```
where POSITION is 0-based, N1 and N2 are the most abundant and the second most abundant nucleotides at this position, respectively.

-  compute_pca_clusters.py: takes the file with R-ratios of GSVs and compute clusters of subjects using PCA:
```
python compute_pca_clusters.py gsv_r_ratios.txt output_dir
```

### Python dependencies:
- scipy, matplotlib, numpy (usually a part of the standard python installation with conda)
- kneed (conda install -c conda-forge kneed)
- pandas (conda install -c anaconda pandas)
- seaborn (conda install -c anaconda seaborn)
- sklearn (conda install -c anaconda scikit-learn)

## Reference:
Safonova Y, Shin SB, Kramer L, Reecy J, Watson CT, Smith TPL, Pevzner PA, Revealing how variations in antibody repertoires correlate with vaccine responses, 2021 [[preprint](https://www.biorxiv.org/content/10.1101/2021.08.06.454618v1)].
