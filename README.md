# parse_mykrobe_predict.py

This script parses Mykrobe predict results for *Shigella sonnei*. 

Mykrobe v0.9.0+ can identify input genomes as _S. sonnei_, assign those identified as _S. sonnei_ to hierarchical genotypes based on detection of single nucleotide variants (SNVs; defined in the file alleles.txt), and report known mutations in the quinolone-resistance determining region (QRDR) of genes *gyrA* (S83L, D87G, D87Y) and *parC* (S80I). Details of the genotyping scheme are available in the preprint [Hawkey et al, 2020](https://www.biorxiv.org/content/10.1101/2020.10.29.360040v1).

This script can be used to parse the resulting JSON files output by Mykrobe (one per genome), and tabulate the results in a single tab-delimited file.

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

| genome     | species   | spp_confidence | final_genotype | name              | num_qrdr | parC_S80I | gyrA_S83L | gyrA_S83A | gyrA_D87G | gyrA_D87N | gyrA_D87Y | num_good_nodes | all_genotype_calls           |
|------------|-----------|----------------|----------------|-------------------|----------|-----------|-----------|-----------|-----------|-----------|-----------|----------------|------------------------------|
| SRR3441855 | S. sonnei | 95.276         | 3.7.18         | Global III        | 0        | 0         | 0         | 0         | 0         | 0         | 0         | 3              | lineage3.7.18                |
| SRR3441874 | S. sonnei | 95.311         | 3.6.1.1        | CipR              | 3        | 1         | 1         | 0         | 1         | 0         | 0         | 4              | lineage3.7.20;lineage3.6.1.1 |
| SRR3446557 | S. sonnei | 95.205         | 2.7.1          | -                 | 1        | 0         | 1         | 0         | 0         | 0         | 0         | 3              | lineage2.7.1                 |
| SRR3452053 | S. sonnei | 95.291         | 2.12           | Latin America IIb | 0        | 0         | 0         | 0         | 0         | 0         | 0         | 2              | lineage2.12                  |

Explanation of columns in the output:
* _genome_: sample ID
* _species_: species call from Mykrobe (_S. sonnei_ or unknown; determined by matching to the ST152 complex from the Achtman _E. coli_ 7-locus MLST scheme)
* _spp_confidence_: percent match to species markers (>90% indicates _S. sonnei_)
* _final_genotype_: final genotype call from Mykrobe, using the scheme of [Hawkey et al, 2020](https://www.biorxiv.org/content/10.1101/2020.10.29.360040v1).
* _name_: human readable alias for genotype, where available (e.g. Global III defined in Holt et al, Nat Genet 2012 corresponds to clade 3.7 in this scheme, so 'Global III' appears in the 'name' column for all genotypes belonging to clade 3.7; see [Hawkey et al, 2020](https://www.biorxiv.org/content/10.1101/2020.10.29.360040v1) for details).
* _num_qrdr_: Total number of mutations detected in the quinolone-resistance determining regions (QRDR) of genes _gyrA_ and _parC_
* _parC_S80I_, _gyrA_S83L_, _gyrA_S83A_, _gyrA_D87G_, _gyrA_D87N_, _gyrA_D87Y_: calls for each individual QRDR mutation. 0 indicates mutation is absent, 1 indicates mutation is present
* _num_good_nodes_: Total number of nodes in the genotyping hierarchy that support the final genotype. For example, for genotype 3.7.18, a num_good_nodes value of 3 indicates a confident call for the single nucleotide variant (SNV) defining lineage 3, a confident call for the SNV defining clade 3.7, and a confident call for the SNV defining subclade 3.7.18.
* _all_genotype_calls_: All genotype calls made by Mykrobe. For more details on each individual call, consult the original Mykrobe json output (under the section labelled "lineage")

#### Unexpected results
As recombination is extremely rare in _S. sonnei_, it is unlikely that DNA isolated from a pure culture would have high quality SNVs that are compatible with multiple genotypes, or that are incompatible with the defined hierarchical scheme of [Hawkey et al, 2020](https://www.biorxiv.org/content/10.1101/2020.10.29.360040v1). Therefore if you see incompatible genotypes reported in the last column, it is worth investigating further whether you have contaminated sequence data or a genuine recombination. A good place to check this is in the JSON output from Mykrobe, under the section labelled "lineage". Each level in the genotyping hierarchy will have a call of either 0 (no kmers found), 1 (confident call) or 0.5 (mixed call). If you see many 0.5 calls, this may indicate a mixed or recombinant strain. An explanation of the JSON output can be found in the [Mykrobe documentation here](https://github.com/Mykrobe-tools/mykrobe/wiki/AMR-prediction-output#json-file).
