# TreeGrafter

This repository contains the code for InterPro's implementation of TreeGrafter (1).

Unlike the [original implementation](https://github.com/pantherdb/TreeGrafter), we use [EPA-ng](https://github.com/Pbdas/epa-ng) (2) instead of [RAxML](https://github.com/stamatak/standard-RAxML) to graft the sequence to the annotated family tree.

## Getting started

Download the PANTHER data, and pre-process them:

```bash
wget http://data.pantherdb.org/ftp/downloads/TreeGrafter/PANTHER16.0_data.tar.gz
tar -zxvf PANTHER16.0_data.tar.gz
python treegrafter/preprocess.py /hps/nobackup/production/interpro/mblum/panther/PANTHER16.0_data/
```

Run hmmsearch (3) on your input sequences:

```
hmmpress PANTHER16.0_data/famhmm/binHmm
hmmsearch PANTHER16.0_data/famhmm/binHmm query.fasta > hits.out
```

Then, run TreeGrafter. TreeGrafter takes three arguments as input (the query sequence file, the hmmsearch output file, and the directory of PANTHER data) and produces a TSV file.
  
```
python treegrafter/treegrafter.py query.fasta hits.out PANTHER16.0_data > predictions.out
```

The options are:

| Option   | Description                                      |
| -------- | ------------------------------------------------ |
| -e       | e-value cutoff (default: disabled)               |
| -o       | output file (instead of the standard output)     |
| --epa-ng | path to the EPA-ng binary (if not in PATH)       |
| -t       | number of threads for EPA-ng to use (default: 1) |
| -T       | path where a temporary directory is created      |
| --keep   | keep temporary directory (default: disabled)     |

The columns of the output TSV are:

| Col | Type    | Description                                      |
| --- | ------- | ------------------------------------------------ |
| 1   | string  | Query ID |
| 2   | string  | Best matched PANTHER family |
| 3   | string  | Predicted PANTHER subfamily |
| 4   | string  | Node of the reference tree where the sequence was grafted onto |
| 5   | float   | Sequence bit score |
| 6   | float   | Sequence E-Value |
| 7   | float   | Domain bit score |
| 8   | float   | Domain E-Value |
| 9   | integer | Start of local alignment (respect to the query profile) |
| 10  | integer | End of local alignment start (respect to the query profile)  |
| 11  | integer | Start of local alignment (respect to the target sequence) |
| 12  | integer | End of local alignment start (respect to the target sequence)  |
| 13  | integer | Start of the envelope of the domain's location (on the target sequence) |
| 14  | integer | End of the envelope of the domain's location (on the target sequence) |
| 15  | string  | PAINT annotations |


### References

1. Haiming Tang, Robert D Finn, Paul D Thomas, TreeGrafter: phylogenetic tree-based annotation of proteins with Gene Ontology terms and other annotations, _Bioinformatics_, Volume 35, Issue 3, 01 February 2019, Pages 518–520, https://doi.org/10.1093/bioinformatics/bty625
2. Pierre Barbera, Alexey M Kozlov, Lucas Czech, Benoit Morel, Diego Darriba, Tomáš Flouri, Alexandros Stamatakis, EPA-ng: Massively Parallel Evolutionary Placement of Genetic Sequences, _Systematic Biology_, Volume 68, Issue 2, March 2019, Pages 365–369, https://doi.org/10.1093/sysbio/syy054
3. http://hmmer.org/
