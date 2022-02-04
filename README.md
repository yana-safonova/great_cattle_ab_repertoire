# The great cattle antibody repertoire
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
- IGHV gene hit, IGHJ gene hit,
- Nucleotide and amino acid sequences of CDR1, CDR2, CDR3,
- Whether or not CDR3 is ultralong (i.e., the length of CDR3 exceeds 150 nt),
- ID of the clone. Clones are computed for all sequences from a given time point. A clone is defined as a group of sequences sharing V hit and J hit and having identical CDR3s.
- Alignment of the variable gene segment against the closest germline IGHV gene.

## Reference:
Safonova Y, Shin SB, Kramer L, Reecy J, Watson CT, Smith TPL, Pevzner PA, Revealing how variations in antibody repertoires correlate with vaccine responses, 2021 [[preprint](https://www.biorxiv.org/content/10.1101/2021.08.06.454618v1)].
