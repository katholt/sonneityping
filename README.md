# parse_mykrobe_predict.py

This script parses Mykrobe predict results for *Shigella sonnei*. 

Mykrobe v0.9.0+ can identify input genomes as _S. sonnei_, assign those identified as _S. sonnei_ to hierarchical genotypes based on detection of single nucleotide variants (SNVs; defined in the file alleles.txt), and report known mutations in the quinolone-resistance determining region (QRDR) of genes *gyrA* (S83L, D87G, D87Y) and *parC* (S80I). Details of the genotyping scheme are available in the preprint [Hawkey et al, 2020](https://www.biorxiv.org/content/10.1101/2020.10.29.360040v1).

This script can be used to parse the resulting JSON files output by Mykrobe (one per genome), and tabulate the results in a single tab-delimited file [(example below)](#example-output).

## Dependencies for parser script
* Python 3.7+
* pandas

## Usage

### Install Mykrobe
First, install Mykrobe (v0.9.0+) as per the instructions on the [Mykrobe github](https://github.com/Mykrobe-tools/mykrobe).

Once Mykrobe is installed, make sure you run the following two commands to ensure you have the most up-to-date panels for genotyping:
```
mykrobe panels update_metadata
mykrobe panels update_species_all
```

### Run mykrobe predict on each genome (Illumina reads, nanopore reads, etc)

Example command (using Illumina reads, SRR3441855):
```
mykrobe predict SRR3441855 sonnei --format json --out SRR3441855.json --seq SRR3441855_1.fastq.gz SRR3441855_2.fastq.gz
```

* For Oxford Nanopore reads, add the flag `--ont` to your command.

* For full details on all Mykrobe options, please see [the Mykrobe documentation](https://github.com/Mykrobe-tools/mykrobe).

### Parse Mykrobe output

**Input**
* JSON files output from `mykrobe predict` using the `sonnei` panel (`--jsons`)

**Output**
* tab-delimited file with one row per genome detailing genotype and QRDR mutations (`--prefix`)

**Example command**
```
python parse_mykrobe_predict.py --jsons mykrobe_results/*.json --prefix results_mykrobe_parsed
```

## Example output
The output table will be named *prefix*_predictResults.tsv, and will be in tab-delimited format:

| genome     | species   | final genotype | name      | confidence        | num QRDR | parC_S80I | gyrA_S83L | gyrA_S83A | gyrA_D87G | gyrA_D87N | gyrA_D87Y | lowest support for genotype marker | poorly supported markers      | max support for additional markers | additional markers         | node support                                                                                                                    |
|------------|-----------|----------------|-----------|-------------------|----------|-----------|-----------|-----------|-----------|-----------|-----------|------------------------------------|-------------------------------|------------------------------------|----------------------------|---------------------------------------------------------------------------------------------------------------------------------|
| sampleA | S.sonnei | 3.6.1.1        | CipR      | strong            | 3        | 1         | 1         | 0         | 1         | 0         | 0         |                                    |                               |                                    |                            | lineage3 (1; 97/0); lineage3.6 (1; 120/0); lineage3.6.1 (1; 91/0);   lineage3.6.1.1 (1; 96/0)                                   |
| sampleB | S.sonnei | 3.6.1.1        | CipR      | strong            | 3        | 1         | 1         | 0         | 1         | 0         | 0         |                                    |                               | 0.009                              | lineage3.7.20 (0.5; 1/105) | lineage3 (1; 95/0); lineage3.6 (1; 112/0); lineage3.6.1 (1; 89/0);   lineage3.6.1.1 (1; 111/1)                                  |
| sampleC | S.sonnei | 3.6.1.1.2      | CipR.MSM5 | moderate          | 3        | 1         | 1         | 0         | 1         | 0         | 0         | 0.792                              | lineage3.6.1.1.2 (0.5; 84/22) |                                    |                            | lineage3 (1; 113/0); lineage3.6 (1; 138/0); lineage3.6.1 (1; 100/0);   lineage3.6.1.1 (1; 131/0); lineage3.6.1.1.2 (0.5; 84/22) |

Explanation of columns in the output:
* **genome**: sample ID
* **species**: species call from Mykrobe (_S. sonnei_ or unknown; determined by matching to the ST152 complex from the Achtman _E. coli_ 7-locus MLST scheme)
* **final genotype**: final genotype call from Mykrobe, using the scheme of [Hawkey et al, 2020](https://www.biorxiv.org/content/10.1101/2020.10.29.360040v1).
* **name**: human readable alias for genotype, where available (e.g. Global III defined in Holt et al, Nat Genet 2012 corresponds to clade 3.7 in this scheme, so 'Global III' appears in the 'name' column for all genotypes belonging to clade 3.7; see [Hawkey et al, 2020](https://www.biorxiv.org/content/10.1101/2020.10.29.360040v1) for details).
* **confidence**: measure of confidence in the final genotype call, summarising read support across all levels in the hierarchy (lineage, clade, subclade, etc)
  * _strong_ - high quality calls (quality of '1' reported by Mykrobe) for ALL levels;
  * _moderate_ - reduced confidence for ONE node (Mykrobe quality of '0.5', but still with >50% read support for the marker allele), high quality calls for ALL OTHERS;
  * _weak_ - low quality for one or more nodes (Mykrobe quality of '0.5' and <50% read support OR Mykrobe quality of '0').
* **num QRDR**: Total number of mutations detected in the quinolone-resistance determining regions (QRDR) of genes _gyrA_ and _parC_
* **parC_S80I, gyrA_S83L, gyrA_S83A, gyrA_D87G, gyrA_D87N, gyrA_D87Y**: calls for each individual QRDR mutation. 0 indicates mutation is absent, 1 indicates mutation is present.
* **lowest support for genotype marker**: For any markers in the final genotype call that do not have a Mykrobe quality of '1', this column reports the percentage of reads supporting the marker allele at the most poorly supported marker (details of all such markers appear in the 'poorly supported markers column').
  * In the example table above, sample sampleC has a low quality call for the final marker in the hierarchy, 3.6.1.1.2, where 79.2% of reads matching this marker belong to the alternate allele.
* **poorly supported markers**: Lists any markers in the final genotype call that do not have Mykrobe quality of '1'. Markers are separated by ';', values in brackets represent the quality call from Mykrobe, followed by the read depth at the alternate / reference alleles. The lowest read support amongst these markers is reported in the previous column.
  * In the example table above, lineage3.6.1.1.2 (0.5, 84/22) is a poorly supported marker for sampleC. The lineage3.6.1.1.2 marker has a Mykrobe quality of 0.5, with a read depth at the marker locus of 84 for the marker allele vs 22 for the reference (wild-type) allele. If no read depth is indicated and only (0) is printed after the marker, this indicates that no kmers matching the marker locus were found in the genome.
* **max support for additional markers**: For any markers detected that are incongruent with the finall genotype call, this column reports the percentage of reads supporting the marker allele at the best supported additional marker.
  * In the example above, sample sampleB is being called as genotype 3.6.1.1 - compatible markers supporting this are 3, 3.6, 3.6.1, 3.6.1.1. However Mykrobe also detected an additional marker (3.7.20, details reported in the 'additional markers' column), with 0.9% read support.
* **additional markers**: Lists any markers that are incongruent with the final genotype call. Markers are separated by ';', and the format is identical to column _poorly supported markers_. The highest read support for any such marker is reported in the previous column.
  * For sampleB, Mykrobe has detected a lineage3.7.20 marker, with quality 0.5. The number of reads matching the marker allele is 1, vs 105 reads matching the reference allele, hence reported as "lineage3.7.20 (0.5, 1/105)".
* **node support**: A list of all markers in the final genotype call with their Mykrobe quality calls (1, 0.5, or 0) and the read depths at the marker allele / reference allele (as per _poorly supported markers_ and _additional markers_). 

**Please note that there is no marker for lineage5.1, as all lineage 5 strains cluster into the same clade (5.1), so this marker is never reported in the output.**

### Unexpected results
As recombination is extremely rare in _S. sonnei_, it is unlikely that DNA isolated from a pure culture would have high quality SNVs that are compatible with multiple genotypes, or that are incompatible with the defined hierarchical scheme of [Hawkey et al, 2020](https://www.biorxiv.org/content/10.1101/2020.10.29.360040v1). Therefore if you see anything reported in the "additional markers" column, it is worth investigating further whether you have contaminated sequence data (which would produce low-quality calls with low read support at additional markers) or a genuine recombination (with would produce high-quality calls with high read support at additional markers). In such cases it may also be illuminating to investigate the Mykrobe output file for further information, details of the JSON output can be found in the [Mykrobe documentation here](https://github.com/Mykrobe-tools/mykrobe/wiki/AMR-prediction-output#json-file).
