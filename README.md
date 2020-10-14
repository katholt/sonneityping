# parse_mykrobe_predict.py

This script parsers Mykrobe predict output for *Shigella sonnei* results, which will assign genomes to genotypes and detect known mutations in the quinolone-resistance determining region (QRDR) of genes *gyrA* (S83L, D87G, D87Y) and *parC* (S80I).

In order for this script to work, the output from Mykrobe **must** be in json format.

## Dependencies for parser script
* Python 3.7+
* pandas

## Run Mykrobe predict
First, install Mykrobe (v0.9.0+) as per the instructions on the [Mykrobe github](https://github.com/Mykrobe-tools/mykrobe).

Once Mykrobe is installed, don't forget to run the following two commands to ensure you have the most up-to-date panels for genotyping:
```
mykrobe panels update_metadata
mykrobe panels update_species_all
```

An example command for running Mykrobe on a sonnei genome (in this case using short-reads for sample SRR3441855):
```
mykrobe predict SRR3441855 sonnei --format json --out SRR3441855.json --seq SRR3441855_1.fastq.gz SRR3441855_2.fastq.gz
```

As input, Mykrobe can take either short-reads, fasta sequences, or Oxford Nanopore reads. If using Oxford Nanopore, add the flag `--ont` to your command.

For full details on all Mykrobe options, please see [the Mykrobe documentation](https://github.com/Mykrobe-tools/mykrobe).

Run Mykrobe predict on all your genomes, and then parse the resulting json files with this script (details below)

## Parsing Mykrobe output

**Input**
* json files output from `mykrobe predict` using the `sonnei` panel (`--jsons`)

**Output**
* tab-delimited file with one row per genome detailing genotype and QRDR mutations (`--prefix`)

```
python parse_mykrobe_predict.py --jsons mykrobe_results/*.json --prefix results_mykrobe_parsed
```

### Example output
The output table will be named *prefix*_predictResults.tsv, and will be in tab-delimited format:

| genome     | species   | spp_confidence | final_genotype | name              | num_qrdr | parC_S80I | gyrA_S83L | gyrA_S83A | gyrA_D87G | gyrA_D87N | gyrA_D87Y | num_good_nodes | all_genotype_calls           |
|------------|-----------|----------------|----------------|-------------------|----------|-----------|-----------|-----------|-----------|-----------|-----------|----------------|------------------------------|
| SRR3441855 | S. sonnei | 95.276         | 3.7.18         | Global III        | 0        | 0         | 0         | 0         | 0         | 0         | 0         | 3              | lineage3.7.18                |
| SRR3441874 | S. sonnei | 95.311         | 3.6.1.1        | CipR              | 3        | 1         | 1         | 0         | 1         | 0         | 0         | 4              | lineage3.7.20;lineage3.6.1.1 |
| SRR3446557 | S. sonnei | 95.205         | 2.7.1          | -                 | 1        | 0         | 1         | 0         | 0         | 0         | 0         | 3              | lineage2.7.1                 |
| SRR3452053 | S. sonnei | 95.291         | 2.12           | Latin America IIb | 0        | 0         | 0         | 0         | 0         | 0         | 0         | 2              | lineage2.12                  |

Explanation of columns in the output:
* _genome_: sample ID
* _species_: species call from Mykrobe (determined by matching to the ST152 complex from the _E. coli_ Achtman MLST scheme) - will either be _S. sonnei_ or unknown
* _spp_confidence_: percent match to species markers (>90% indicates _S. sonnei_)
* _final_genotype_: final genotype call from mykrobe
* _name_: human readable alias for genotype, where available (e.g. Global III defined in Holt et al, Nat Genet 2012 corresponds to clade 3.7 in this scheme, so 'Global III' appears in the 'name' column for all genotypes belonging to clade 3.7).
* _num_qrdr_: Total number of QRDR mutations detected
* _parC_S80I_, _gyrA_S83L_, _gyrA_S83A_, _gyrA_D87G_, _gyrA_D87N_, _gyrA_D87Y_: calls for each individual QRDR mutation. 0 indicates mutation is absent, 1 indicates mutation is present
* _num_good_nodes_: Total number of nodes in the genotyping hierarchy that support the final genotype. For example, for genotype 3.7.18, a num_good_nodes value of 3 indicates a confident call for the SNV defining lineage 3, a confident call for the SNV defining clade 3.7, and a confident call for the SNV defining subclade 3.7.18
* _all_genotype_calls_: All genotype calls made by Mykrobe. For more details on each individual call, consult the original Mykrobe json output (under the section labelled "lineage")

#### Unexpected results
As recombination is extremely rare in S. sonnei, it is unlikely that DNA isolated from a single S. sonnei culture would have high quality SNPs that are compatible with multiple genotypes, or incompatible clade/subclade combinations... therefore if you see this in the output, it is worth investigating further whether you have contaminated sequence data or a genuine recombination. A good place to check this is in the json output from Mykrobe, under the section labelled "lineage". Each level in the genotyping hierarchy will have a call of either 0 (no kmers found), 1 (confident call) or 0.5 (mixed call). If you see many 0.5 calls, this may indicate a mixed or recombinant strain.
