# TreeGrafter

This repository contains the code for InterPro's implementation of TreeGrafter (1).

Unlike the [original implementation](https://github.com/pantherdb/TreeGrafter), we use [EPA-ng](https://github.com/Pbdas/epa-ng) (2) instead of [RAxML](https://github.com/stamatak/standard-RAxML) (3) to graft the sequence to the annotated family tree.

## Getting started

TreeGrafter now uses two separate data directories:

- **Library directory** — the PANTHER library (trees, alignments, HMMs). Updated yearly.
- **Annotation directory** — per-family JSON files derived from PAINT annotations. Updated monthly. (More info [here](https://github.com/pantherdb/fullgo_paint_update/issues/77))

### 1. Download PANTHER library and PAINT annotations

```bash
$ wget http://data.pantherdb.org/ftp/downloads/TreeGrafter/PANTHER19.0_data.tar.gz
$ tar -zxvf PANTHER19.0_data.tar.gz
$ wget https://data.pantherdb.org/ftp/downloads/paint/19.0/2026-01-05/PAINT_TreeGrafter_Annotations_TOTAL.txt.gz
$ gunzip PAINT_TreeGrafter_Annotations_TOTAL.txt.gz
```

### 2. Prepare annotations

Convert the PAINT annotation file into per-family JSONs. This is only required once per annotation release:

```bash
$ python treegrafter.py prepare PAINT_TreeGrafter_Annotations_TOTAL.txt annotations/
```

### 3. Run hmmsearch

Run hmmsearch (4) on your input sequences:

```
$ hmmsearch PANTHER19.0_data/famhmm/binHmm query.fasta > hits.out
```

### 4. Run TreeGrafter

```
$ python treegrafter.py run query.fasta hits.out -d PANTHER19.0_data -a annotations/ > predictions.tsv
```

### Options

When running `treegrafter.py run`, options are:

| Option   | Description                                      |
| -------- | ------------------------------------------------ |
| -d       | **required** — PANTHER library directory         |
| -a       | **required** — annotation directory (from `prepare`) |
| -e       | e-value cutoff (default: disabled)               |
| -o       | output file (instead of the standard output)     |
| --epa-ng | path to the EPA-ng binary (if not in PATH)       |
| -t       | number of threads for EPA-ng to use (default: 1) |
| -T       | path where a temporary directory is created      |
| --keep   | keep temporary directory (default: disabled)     |
| --print-go | include GO terms and protein class columns in output |

### Output format

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
| 14  | string  | GO terms (only with `--print-go`) |
| 15  | string  | PANTHER protein class (only with `--print-go`) |

## Docker

TreeGrafter is available as a Docker image. PANTHER data and annotations need to be provided to the container with bind mounts. Assuming both directories are under your current working directory, use `-v $(pwd):/mnt`.

To prepare annotations:

```bash
$ docker run --rm -v "$(pwd)":/mnt interpro/treegrafter prepare /mnt/PAINT_Annotations_TOTAL.txt /mnt/annotations
```

To search your sequences:

```bash
$ docker run --rm -v "$(pwd)":/mnt interpro/treegrafter search /mnt/query.fasta /mnt/PANTHER19.0_data /mnt/annotations /mnt/predictions.tsv
```

## References

1. Haiming Tang, Robert D Finn, Paul D Thomas, TreeGrafter: phylogenetic tree-based annotation of proteins with Gene Ontology terms and other annotations, _Bioinformatics_, Volume 35, Issue 3, February 2019, Pages 518–520, https://doi.org/10.1093/bioinformatics/bty625
2. Pierre Barbera, Alexey M Kozlov, Lucas Czech, Benoit Morel, Diego Darriba, Tomáš Flouri, Alexandros Stamatakis, EPA-ng: Massively Parallel Evolutionary Placement of Genetic Sequences, _Systematic Biology_, Volume 68, Issue 2, March 2019, Pages 365–369, https://doi.org/10.1093/sysbio/syy054
3. Alexandros Stamatakis, RAxML version 8: a tool for phylogenetic analysis and post-analysis of large phylogenies, _Bioinformatics_, Volume 30, Issue 9, May 2014, Pages 1312–1313, https://doi.org/10.1093/bioinformatics/btu033
4. http://hmmer.org/
