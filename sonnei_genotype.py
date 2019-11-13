#!/usr/bin/env python
#
# Input BAM (recommended) or VCF (if highly trusted SNP data) relative to Shigella sonnei 53G (NC_016822) and assign S. sonnei genotype codes and detect QRDR mutations.
#
# Authors - Kat Holt (drkatholt@gmail.com), Kalani Paranagama
#
# Documentation - https://github.com/katholt/sonneityping
#
# Dependencies:
#	 SAMtools (v1.2) and bcftools (v1.2) are required to genotype from BAMs.
#
# Last modified - Sep 27, 2019
#

from argparse import (ArgumentParser, FileType)
import os, sys, re, collections, operator
import gzip
import logging
import time
import datetime
from subprocess import call, check_output, CalledProcessError, STDOUT


def parse_args():
	"Parse the input arguments, use '-h' for help"
	parser = ArgumentParser(description='VCF to Shigella sonnei genotypes')
	parser.add_argument('--mode', required=True,
						help='Mode to run in based on input files (vcf, bam, or vcf_parsnp)')
	parser.add_argument(
		'--vcf', nargs='+', type=str, required=False,
		help='VCF file(s) to genotype (Mapping MUST have been done using 53G as a reference sequence)')
	parser.add_argument('--bam', nargs='+', type=str, required=False,
						help='BAM file(s) to genotype (Mapping MUST have been done using 53G as a reference sequence')
	parser.add_argument('--ref_id', type=str, required=False,
						help='Name of the reference in the VCF file (#CHROM column)')
	parser.add_argument('--phred', type=int, required=False, default=20,
						help='Minimum phred quality to count a variant call vs 53G as a true SNP (default 20)')
	parser.add_argument('--min_prop', type=float, required=False, default=0.1,
						help='Minimum proportion of reads required to call a SNP (default 0.1)')
	parser.add_argument('--ref', type=str, required=False,
						help='Reference sequence in fasta format. Required if bam files provided.')
	parser.add_argument('--output', type=str, required=False, default='genotypes.txt',
						help='Location and name for output file.')
	parser.add_argument('--samtools_location', type=str, required=False, default='',
						help='Location of folder containing samtools installation if not standard/in path.')
	parser.add_argument('--bcftools_location', type=str, required=False, default='',
						help='Location of folder containing bcftools installation if not standard/in path.')
	parser.add_argument('--allele_table', type=str, required=False, default='',
						help='Location of allele table.')
	return parser.parse_args()



### QRDR SNP definitions

qrdr_loci = [2503749,2503737, 2503738, 3476068, 2503750, 2503738]
qrdr_snp_alleles = ['A', 'C', 'A', 'A', 'C', 'T']
qrdr_groups = ['gyrA-S83L', 'gyrA-D87G', 'gyrA-D87Y', 'parC-S80I', 'gyrA-S83A', 'gyrA-D87N']


### Genotype SNP definitions

loci = []
snp_alleles = []
groups = []
group_names = {} #key=group, value=name


# check if this SNP defines a QRDR group
def checkQRDRSNP(vcf_line_split, this_qrdr_groups, qrdr_proportions, args):
	qrdr_snp = int(vcf_line_split[1])
	if qrdr_snp in qrdr_loci:
		i = qrdr_loci.index(qrdr_snp)
		if float(vcf_line_split[5]) > args.phred:
			m = re.search("DP4=(\d+),(\d+),(\d+),(\d+)", vcf_line_split[7])
			if m != None:
				alt_read_count = int(m.group(3)) + int(m.group(4))
				total_read_count = alt_read_count + int(m.group(1)) + int(m.group(2))
				if float(total_read_count) == 0:
					qrdr_snp_proportion = float(-1)
				else:
					qrdr_snp_proportion = float(alt_read_count) / total_read_count
			else:
				if vcf_line_split[4] != '.':  # if the ALT is not '.' i.e. if the alt is not same as ref
					try:
						ad = vcf_line_split[9].split(':')[1].split(',')	 # get the AD ratio
						alt_read_count = int(ad[1])
						total_read_count = int(ad[0]) + alt_read_count
						qrdr_snp_proportion = float(alt_read_count) / total_read_count
					except IndexError:
						qrdr_snp_proportion = float(-1)

				else:
					qrdr_snp_proportion = float(
						-1)	 # set unknowns to negative so that we know this is not a real proportion

			qrdr_snp_allele = vcf_line_split[4]
			for position in xrange(0, len(qrdr_loci)):
				if (qrdr_snp == qrdr_loci[position]) and (qrdr_snp_allele == qrdr_snp_alleles[position]) and (
						qrdr_snp_proportion > args.min_prop):
					this_qrdr_groups.append(qrdr_groups[position])	# Add QRDR SNP

	return (this_qrdr_groups)


# check if this SNP defines a group
def checkSNP(vcf_line_split, this_groups, proportions, args):
	snp = int(vcf_line_split[1])
	if snp in loci:
		i = loci.index(snp)

		if float(vcf_line_split[5]) > args.phred:
			m = re.search("DP4=(\d+),(\d+),(\d+),(\d+)", vcf_line_split[7])
			if m != None:
				alt_read_count = int(m.group(3)) + int(m.group(4))
				total_read_count = alt_read_count + int(m.group(1)) + int(m.group(2))
				if float(total_read_count) == 0:
					snp_proportion = float(-1)
				else:
					snp_proportion = float(alt_read_count) / total_read_count
			else:
				if vcf_line_split[4] != '.':  # if the ALT is not '.' i.e. if the alt is not same as ref
					try:
						ad = vcf_line_split[9].split(':')[1].split(',')	 # get the AD ratio
						alt_read_count = int(ad[1])
						total_read_count = int(ad[0]) + alt_read_count
						snp_proportion = float(alt_read_count) / total_read_count
					except IndexError:
						snp_proportion = float(-1)

				else:
					snp_proportion = float(-1)	# set unknowns to negative so that we know this is not a real proportion
			if snp_proportion > args.min_prop:
				this_allele = vcf_line_split[4]
				if this_allele == snp_alleles[i]:
					this_groups.append(groups[i])  # retrieve the group that this SNP defines
					proportions[groups[i]] = snp_proportion
	return (this_groups, proportions)


def checkSNPmulti(vcf_line_split, this_groups, args):
	snp = int(vcf_line_split[1])
	if snp in loci:
		i = loci.index(snp)
		strain = 0
		for gt in vcf_line_split[10:]:
			if (int(gt) == 1) and (vcf_line_split[4] == snp_alleles[i]):
				if strain in this_groups:
					this_groups[strain].append(groups[i])  # retrieve the group that this SNP defines
				else:
					this_groups[strain] = [groups[i]]
			strain += 1
	return this_groups


# sort groups into the three levels (primary, clade, subclade)
def parseGeno(this_groups, proportions, group_names):

	subclades = []
	clades = []
	primary = []
	for group in this_groups:
		level = len(group.split("."))
		if level >= 3:
			subclades.append(group)
		elif level == 2:
			clades.append(group)
		elif level == 1:
			primary.append(group)

#state clade for singleton clades
	if ('1.1.1' in subclades):
		clades.append('1.1')
	if ('1.2.1') in subclades:
		clades.append('1.2')
	if ('1.3.1') in subclades:
		clades.append('1.3')
	if ('1.4.1' in subclades):
		clades.append('1.4')
	if ('2.2.1' in subclades):
		clades.append('2.2')
	if ('2.3.1' in subclades):
		clades.append('2.3')
	if ('3.8.1' in subclades):
		clades.append('3.8')
	if ('3.9.1' in subclades):
		clades.append('3.9')

#state clade for single-subclade clades
	if ('3.1.1' in subclades):
		clades.append('3.1')
	if ('3.2.1' in subclades):
		clades.append('3.2')
	if ('3.5.1' in subclades):
		clades.append('3.5')

#state clade for lineage 5
	if ('5.1.1' in subclades):
		clades.append('5.1')
	if ('5.1.2' in subclades):
		clades.append('5.1')
	if ('5.1.3' in subclades):
		clades.append('5.1')
	if ('5.1.4' in subclades):
		clades.append('5.1')
	if ('5.1.5' in subclades):
		clades.append('5.1')
	if ('5.1.6' in subclades):
		clades.append('5.1')

#fix 3.3/3.4 nesting
	if ('3.4.1' in subclades) and ('3.4' in clades) and ('3.3' in clades):
		clades.remove('3.3')
	if ('3.4.2' in subclades) and ('3.4' in clades) and ('3.3' in clades):
		clades.remove('3.3')
	if ('3.4.3' in subclades) and ('3.4' in clades) and ('3.3' in clades):
		clades.remove('3.3')
	if ('3.4.4' in subclades) and ('3.4' in clades) and ('3.3' in clades):
		clades.remove('3.3')
	if ('3.4.5' in subclades) and ('3.4' in clades) and ('3.3' in clades):
		clades.remove('3.3')
	if ('3.4.6' in subclades) and ('3.4' in clades) and ('3.3' in clades):
		clades.remove('3.3')

#3.7.3/3.7.2 nesting
	if ('3.7.3' in subclades) and ('3.7.2' in subclades):
		subclades.remove('3.7.3')


#lineage 4 - state lineage and clade
	if ('4.1.1' in subclades):
		clades.append('4.1')


#fix 3.6.1/3.6.1.1/3.6.1.1.1/3.6.1.1.2 nesting
	if ('3.6.1' in subclades) and ('3.6.1.1' in subclades) and ('3.6.1.1.1' in subclades):
		subclades.remove('3.6.1')
		subclades.remove('3.6.1.1')
	if ('3.6.1' in subclades) and ('3.6.1.1' in subclades) and ('3.6.1.1.2' in subclades):
		subclades.remove('3.6.1')
		subclades.remove('3.6.1.1')


#fix 3.6.1.1.3/3.6.1.1.3.1 nesting
	if ('3.6.1' in subclades) and ('3.6.1.1' in subclades) and ('3.6.1.1.3' in subclades) and ('3.6.1.1.3.1' in subclades):
		subclades.remove('3.6.1')
		subclades.remove('3.6.1.1')
		subclades.remove('3.6.1.1.3')
	if ('3.6.1' in subclades) and ('3.6.1.1' in subclades) and ('3.6.1.1.3' in subclades):
		subclades.remove('3.6.1')
		subclades.remove('3.6.1.1')

#fix 3.6.1/3.6.1.1 nesting
	if ('3.6.1'in subclades) and ('3.6.1.1'in subclades):
		subclades.remove('3.6.1')

#fix 3.7.29/3.7.29.1/3.7.29.1.1/3.7.29.1.1.1/3.7.29.1.1.2 nesting
	if ('3.7.29' in subclades) and ('3.7.29.1' in subclades) and ('3.7.29.1.1' in subclades) and ('3.7.29.1.1.1' in subclades):
		subclades.remove('3.7.29')
		subclades.remove('3.7.29.1')
		subclades.remove('3.7.29.1.1')
	if ('3.7.29' in subclades) and ('3.7.29.1' in subclades) and ('3.7.29.1.1' in subclades) and ('3.7.29.1.1.2' in subclades):
		subclades.remove('3.7.29')
		subclades.remove('3.7.29.1')
		subclades.remove('3.7.29.1.1')
	if ('3.7.29' in subclades) and ('3.7.29.1' in subclades) and ('3.7.29.1.1' in subclades):
		subclades.remove('3.7.29')
		subclades.remove('3.7.29.1')


#fix 3.7.29/3.7.29.1/3.7.29.1.2/3.7.29.1.2.1 nesting
	if ('3.7.29' in subclades) and ('3.7.29.1' in subclades) and ('3.7.29.1.2' in subclades) and ('3.7.29.1.2.1' in subclades):
		subclades.remove('3.7.29')
		subclades.remove('3.7.29.1')
		subclades.remove('3.7.29.1.2')
	if ('3.7.29' in subclades) and ('3.7.29.1' in subclades) and ('3.7.29.1.2' in subclades):
		subclades.remove('3.7.29')
		subclades.remove('3.7.29.1')


#fix 3.7.29/3.7.29.1/3.7.29.1.3 nesting
	if ('3.7.29' in subclades) and ('3.7.29.1' in subclades) and ('3.7.29.1.3' in subclades):
		subclades.remove('3.7.29')
		subclades.remove('3.7.29.1')

#fix 3.7.29/3.7.29.1/3.7.29.1.4/3.7.29.1.4.1 nesting
	if ('3.7.29' in subclades) and ('3.7.29.1' in subclades) and ('3.7.29.1.4' in subclades) and ('3.7.29.1.4.1' in subclades):
		subclades.remove('3.7.29')
		subclades.remove('3.7.29.1')
		subclades.remove('3.7.29.1.4')
	if ('3.7.29' in subclades) and ('3.7.29.1' in subclades) and ('3.7.29.1.4' in subclades):
		subclades.remove('3.7.29')
		subclades.remove('3.7.29.1')


#fix 3.7.29/3.7.29.1 nesting
	if ('3.7.29' in subclades) and ('3.7.29.1' in subclades):
		subclades.remove('3.7.29')

#fix 3.7.30.4/3.7.30.4.1 nesting
	if ('3.7.30.4' in subclades) and ('3.7.30.4.1' in subclades):
		subclades.remove('3.7.30.4')


### KH - this is where we need to deal with the fact that the reference is 2.8.2
###		so strains of genotype 2.8.2 should have no SNPs in their VCF relative to
###		the reference which is also 2.8.2


# fix primary clades relative to 53G = 2.8.2, ie has clade1, clade2, clade3 SNPs
	if '2' in primary:
		primary.remove('2')  # anything NOT in 2 will have this SNP
	else:
		if len(primary) == 0:
			primary.append('2')  # anything with no clade, and 2 SNP not called, belongs in primary 2 with 53G
		

# fix clade relative to 53G:
	if '2.8' in clades:
		clades.remove('2.8')  # anything NOT in 2.8 will have this SNP
	else:
		if len(clades) == 0:
			clades.append('2.8')  # anything with no clade, and 2.8 SNP not called, belongs in 2.8 with 53G


# fix subclades relative to 53G:
	if '2.8.2' in subclades:
		subclades.remove('2.8.2')  # anything NOT in 2.8.2 will have this SNP

	# state subclade for 3.6.0 strains
	if (len(subclades) == 0) and ('3.6' in clades):
		subclades.append('3.6.0')

	# state subclade for reference
	if len(subclades) == 0:
		subclades.append('2.8.2')  # anything with no subclade, and 2.8.2 SNP NOT called, belongs in 2.8.2 with 53G



#store final genotype, to the lowest level available
	if len(primary) == 0:
		primary.append(group.split(".")[0])
		
	final_geno = primary[0]
	if len(clades) > 0:
		final_geno = ','.join(clades)
	if len(subclades) > 0:
		final_geno = ','.join(subclades)
		#print final_geno
		#print primary
		#print clades

	# add proportion of reads supporting each of these groups
	p_prod = 1

	p_sub = []
	for group in subclades:
		if group in proportions:
			p_sub.append(str(round(proportions[group], 2)))
			p_prod = p_prod * proportions[group]

	p_cl = []
	for group in clades:
		if group in proportions:
			p_cl.append(str(round(proportions[group], 2)))
			p_prod = p_prod * proportions[group]

	p_pr = []
	for group in primary:
		if group in proportions:
			p_pr.append(str(round(proportions[group], 2)))
			p_prod = p_prod * proportions[group]

	# final call
	info = final_geno + '\t'
	if 'A' in proportions:
		info += 'A'	 # annotate as 'A' to indicate this comes from assembled data and not reads
	else:
		info += str(round(p_prod, 2))  # indicate proportion of reads supporting this call

	# name
	this_name = '-'
	if final_geno in group_names:
		this_name = group_names[final_geno]
	info += '\t' + this_name

	# level calls
	info += '\t' + ','.join(subclades) + '\t' + ','.join(clades) + '\t' + ','.join(primary)

	# level proportions
	info += '\t' + ','.join(p_sub) + '\t' + ','.join(p_cl) + '\t' + ','.join(p_pr)

	return info


# exception to raise if the command we try to run fails for some reason
class CommandError(Exception):
	pass


def run_command(command, **kwargs):
	'Execute a shell command and check the exit status and any O/S exceptions'
	command_str = ' '.join(command)
	logging.info('Running: {}'.format(command_str))
	try:
		exit_status = call(command, **kwargs)
	except OSError as e:
		message = "Command '{}' failed due to O/S error: {}".format(command_str, str(e))
		raise CommandError({"message": message})
	if exit_status != 0:
		message = "Command '{}' failed with non-zero exit status: {}".format(command_str, exit_status)
		raise CommandError({"message": message})


# main function
def main():
	args = parse_args()

	f = file(args.allele_table,"r")
	for line in f:
		fields = line.strip().split('\t')
		locus = int(fields[0])
		allele = fields[1]
		group = fields[2]
		group_name = fields[3]
		loci.append(locus)
		snp_alleles.append(allele)
		groups.append(group)
		group_names[group]=group_name
	f.close()

	if (((args.mode == 'vcf') and args.vcf and args.ref_id) or (
			(args.mode == 'bam') and args.bam and args.ref and args.ref_id) or (
			(args.mode == 'vcf_parsnp') and args.vcf)):

		# Initialise output file and timestamp
		timestamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y.%m.%d_%H.%M.%S_')
		output_file = open(timestamp + args.output, 'w')

		# GENERATE VCFS (1 per strain) FROM BAMS

		if args.mode == 'bam':

			with open(args.ref, 'r') as fasta_file:	 # create SAMtools compatible fasta file
				sequences = fasta_file.read()
				for sequence in sequences.split('>'):
					if args.ref_id in sequence:
						new_header = '>' + args.ref_id
						replacement_index = sequence.find('\n')
						sequence = new_header + sequence[replacement_index:]
						with open('temp_reference.fasta', 'w') as temp_fasta_file:
							temp_fasta_file.write(sequence)

			vcfFiles = []

			# coordinates in zero-base, half-open for SAMtools compatible bed file
			ordered_loci = list(loci + qrdr_loci)
			sorted(ordered_loci)
			temp_bed_file = open(args.ref_id + '.bed', 'w')	 # create temporary bed file for SAMtools
			for locus in ordered_loci:
				temp_bed_file.write(
					args.ref_id + '\t' + str(locus - 1) + '\t' + str(locus) + '\n')	 # write bed file from matrix
			temp_bed_file.close()  # close bedFile

			run_command(['samtools', 'faidx', 'temp_reference.fasta'])	# index fasta file

			for bam in args.bam:
				print
				'bam files supplied, generating vcf file for ' + bam
				if os.path.exists(bam + '.bai') == False:  # index bam file if indexed bam not provided
					run_command([args.samtools_location + 'samtools', 'index', bam])

				run_command(
					[args.samtools_location + 'samtools', 'mpileup', '-q', str(args.phred), '-ugB', '-f',
					 'temp_reference.fasta',
					 '-l', args.ref_id + '.bed', bam, '-o', bam[:-4] + '.output', '-I'])  # detect SNPs

				run_command(
					[args.bcftools_location + 'bcftools', 'call', '-c', bam[:-4] + '.output', '-o',
					 bam[:-4] + '.vcf'])  # generate vcf files
				run_command(['rm', bam[:-4] + '.output'])

				vcfFiles.append(bam[:-4] + '.vcf')	# supply generated vcf file to script

			run_command(['rm', args.ref_id + '.bed'])  # remove temp files
			run_command(['rm', 'temp_reference.fasta'])
			run_command(['rm', 'temp_reference.fasta.fai'])

			args.vcf = vcfFiles

		# PRINT OUTPUT HEADER

		if args.mode == 'bam':
			output_file.write('\t'.join(
				['File', 'Genotype', 'Genotype_support', 'Genotype_name', 'Subclade', 'Clade', 'Lineage', 'Support_Subclade',
				 'Support_Clade', 'Support_Lineage', 'QRDR mutations', 'Number of SNPs called\n']))
		elif args.mode == 'vcf':
			output_file.write('\t'.join(
				['File', 'Genotype', 'Genotype_support', 'Genotype_name', 'Subclade', 'Clade', 'Lineage', 'Support_Subclade',
				 'Support_Clade', 'Support_Lineage', 'QRDR mutations\n']))
		else:
			output_file.write('\t'.join(
				['File', 'Genotype', 'Genotype_support', 'Genotype_name', 'Subclade', 'Clade', 'Lineage', 'Support_Subclade',
				 'Support_Clade', 'Support_Lineage\n']))

		# PARSE MAPPING BASED VCFS (1 per strain)

		if (args.vcf and (args.mode != 'vcf_parsnp')):
			for vcf in args.vcf:
				snp_count = 0
				this_groups = []  # list of groups identified by defining SNPs
				this_qrdr_groups = []  # list of QRDR SNPs found
				proportions = {}  # proportion of reads supporting each defining SNP; key = group, value = proportion
				qrdr_proportions = {}  # proportion of reads supporting each defining SNP; key = group, value = proportion

				# read file
				(file_name, ext) = os.path.splitext(vcf)

				if ext == '.gz':
					f = gzip.open(vcf, 'r')
				else:
					f = open(vcf, 'r')

				any_ref_line = 0

				for line in f:
					if not line.startswith('#'):
						x = line.rstrip().split()
						if not line.strip() == '':
							snp_count = snp_count + 1
						if x[0] == args.ref_id:
							# parse this SNP line
							any_ref_line = 1
							(this_groups, proportions) = checkSNP(x, this_groups, proportions, args)
							this_qrdr_groups = checkQRDRSNP(x, this_qrdr_groups, qrdr_proportions, args)

				f.close()
				

				#"qrdr_groups".join(qrdr_groups)
				if any_ref_line > 0:
					info = parseGeno(this_groups, proportions, group_names)
					if args.bam:
						output_file.write(
							vcf + '\t' + info + '\t' + ','.join(this_qrdr_groups) + '\t' + str(snp_count) + '\n')
					else:
						output_file.write(vcf + '\t' + info + '\t' + ','.join(this_qrdr_groups) + '\n')
				else:
					output_file.write(
						vcf + '\tNo SNPs encountered against expected reference. Wrong reference or no SNP calls?\n')

		# PARSE PARSNP VCF (multiple strains)

		if args.mode == 'vcf_parsnp':

			if not args.ref_id:
				args.ref_id = '1'

			for vcfm in args.vcf:

				# read file
				(file_name, ext) = os.path.splitext(vcfm)

				if ext == '.gz':
					f = gzip.open(vcfm, 'r')
				else:
					f = open(vcfm, 'r')

				any_ref_line = 0

				this_groups = {}  # list of groups identified by defining SNPs, key = strain id (column number)
				strains = []

				for line in f:
					x = line.rstrip().split()
					if x[0] == '#CHROM':
						strains = x[10:]
					if not line.startswith('#'):
						if x[0] == args.ref_id:
							any_ref_line = 1  # parse this SNP line
							this_groups = checkSNPmulti(x, this_groups, args)

				f.close()

				# collate by strain
				if any_ref_line > 0:
					for strain in this_groups:
						info = parseGeno(this_groups[strain], ['A'], group_names)
						output_file.write(strains[strain] + '\t' + info + '\n')
				else:
					output_file.write(strains[strain] + '\tNo SNPs encountered against expected reference. Wrong reference or no SNP calls?\n')

		output_file.close()
	else:
		print
		'Missing or incomplete input parameters, please check these and try again.'


# call main function
if __name__ == '__main__':
	main()
