# TreeGrafter

This repository contains the code for InterPro's implementation of TreeGrafter (1).

Unlike the [original implementation](https://github.com/pantherdb/TreeGrafter), we use [EPA-ng](https://github.com/Pbdas/epa-ng) (2) instead of [RAxML](https://github.com/stamatak/standard-RAxML) (3) to graft the sequence to the annotated family tree.

## Getting started

Download and extract PANTHER data:

```bash
$ wget http://data.pantherdb.org/ftp/downloads/TreeGrafter/PANTHER17.0_data.tar.gz
$ tar -zxvf PANTHER17.0_data.tar.gz
```

Prepare PANTHER annotations. This is only required once:

```bash
$ python treegrafter/treegrafter.py prepare PANTHER17.0_data
```

Run hmmsearch (4) on your input sequences:

```
$ hmmsearch PANTHER17.0_data/famhmm/binHmm query.fasta > hits.out
```

Run TreeGrafter. TreeGrafter takes at least three arguments as input:

1. the query sequence file
2. the hmmsearch output file
3. the directory of prepared PANTHER data
  
```
$ python treegrafter/treegrafter.py run query.fasta hits.out PANTHER17.0_data > predictions.tsv
```

Options are:

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
| 2   | string  | Predicted PANTHER subfamily (if any) or best matched PANTHER family |
| 3   | float   | Sequence bit score |
| 4   | float   | Sequence E-Value |
| 5   | float   | Domain bit score |
| 6   | float   | Domain E-Value |
| 7   | integer | Start of local alignment (respect to the query profile) |
| 8   | integer | End of local alignment start (respect to the query profile)  |
| 9   | integer | Start of local alignment (respect to the target sequence) |
| 10  | integer | End of local alignment start (respect to the target sequence)  |
| 11  | integer | Start of the envelope of the domain's location (on the target sequence) |
| 12  | integer | End of the envelope of the domain's location (on the target sequence) |
| 13  | string  | Node of the reference tree where the sequence was grafted onto |

## Docker

TreeGrafter is available as a Docker image. PANTHER data need to be provided to the container with bind mounts. Assuming the `PANTHER17.0_data` directory is in your current working directory, you can use `-v $(pwd):/mnt` so the PANTHER data will be mounted in `/mnt/PANTHER17.0_data` in the container.

To prepare PANTHER data:

```bash
$ docker run --rm -v "$(pwd)":/mnt interpro/treegrafter prepare /mnt/PANTHER17.0_data
```

To search your sequences:

```bash
$ docker run --rm -v "$(pwd)":/mnt interpro/treegrafter search /mnt/query.fasta /mnt/PANTHER17.0_data /mnt/predictions.tsv
```

## References

1. Haiming Tang, Robert D Finn, Paul D Thomas, TreeGrafter: phylogenetic tree-based annotation of proteins with Gene Ontology terms and other annotations, _Bioinformatics_, Volume 35, Issue 3, February 2019, Pages 518–520, https://doi.org/10.1093/bioinformatics/bty625
2. Pierre Barbera, Alexey M Kozlov, Lucas Czech, Benoit Morel, Diego Darriba, Tomáš Flouri, Alexandros Stamatakis, EPA-ng: Massively Parallel Evolutionary Placement of Genetic Sequences, _Systematic Biology_, Volume 68, Issue 2, March 2019, Pages 365–369, https://doi.org/10.1093/sysbio/syy054
3. Alexandros Stamatakis, RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies, _Bioinformatics_, Volume 30, Issue 9, May 2014, Pages 1312–1313, https://doi.org/10.1093/bioinformatics/btu033
4. http://hmmer.org/
