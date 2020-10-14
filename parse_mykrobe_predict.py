import json
import pandas as pd
from argparse import ArgumentParser

def get_arguments():
    parser = ArgumentParser(description='Parse mykrobe predict JSON files')

    # job submission options
    parser.add_argument('--jsons', required=True, nargs='+', help='JSON files output from mykrobe predict')
    parser.add_argument('--prefix', required=True, help='prefix for output files')

    return parser.parse_args()

def extract_qrdr_info(genome_data, genome_name, spp_call):

    # set up empty dict for final output
    qrdr_out_dict = {}
    # make a list of all possible mutations
    qrdr_possible = ['parC_S80I', 'gyrA_S83L', 'gyrA_S83A', 'gyrA_D87G', 'gyrA_D87N', 'gyrA_D87Y']

    # can only extract qrdr info if sample is sonnei, otherwise set 'unknown' for all calls
    if spp_call == "Shigella_sonnei":
        try:
            qrdr_data = genome_data["susceptibility"]["quinolones"]
        # if nothiing has been called, just set everything to NA
        except KeyError:
            qrdr_out_dict['genome'] = [genome_name]
            for allele in qrdr_possible:
                qrdr_out_dict[allele] = ['NA']
            qrdr_out_dict['num_qrdr'] = ['NA']
            # make table
            qrdr_out_df = pd.DataFrame(qrdr_out_dict, columns=['genome'] + qrdr_possible + ['num_qrdr'])
            return qrdr_out_df
    else:
        qrdr_out_dict['genome'] = [genome_name]
        for allele in qrdr_possible:
            qrdr_out_dict[allele] = ['NA']
        qrdr_out_dict['num_qrdr'] = ['NA']
        # make table
        qrdr_out_df = pd.DataFrame(qrdr_out_dict, columns=['genome'] + qrdr_possible + ['num_qrdr'])
        return qrdr_out_df

    # set up the list of calls in this genome
    qrdr_calls = []
    # only need to parse if the predict is R - if it's S all values will be 0
    if qrdr_data["predict"] == "R":
        calls = qrdr_data["called_by"]
        # loop through each mutation
        for mutation in calls:
            # get the mutation name and add it to our list of calls
            qrdr_calls.append(mutation.split('-')[0])
    # now create our pandas dataframe for this genome
    qrdr_out_dict['genome'] = [genome_name]
    num_qrdr_calls = 0
    for allele in qrdr_possible:
        if allele in qrdr_calls:
            qrdr_out_dict[allele] = [1]
            num_qrdr_calls += 1
        else:
            qrdr_out_dict[allele] = [0]
    # add column with total number of qrdr calls
    qrdr_out_dict['num_qrdr'] = [num_qrdr_calls]
    # make table
    qrdr_out_df = pd.DataFrame(qrdr_out_dict, columns=['genome'] + qrdr_possible + ['num_qrdr'])
    return qrdr_out_df

def inspect_calls(genotype_details):
    # for all the genotypes that mykrobe is calling, inspect the values listed next to the actual call heirarhcy
    genotype_list = list(genotype_details.keys())
    # intialise empty values for best score and corresponding genotype
    best_score = 0
    best_genotype = None
    num_good_nodes = None
    # loop through each genotype
    for genotype in genotype_list:
        # the maximum score we can get is the total depth of the heirarchy
        max_score = genotype_details[genotype]['tree_depth']
        # the actual score is a sum of the values within each levels call
        actual_score = sum(list(genotype_details[genotype]['genotypes'].values()))
        # if actual score is equal to the max score, then this is the best genotype
        # BUT WE NEED TO DEAL WITH INSTANCES WHERE WE MIGHT HAVE A GREAT CALL 3.7.25 say AND A LESS GREAT CALL 3.6.1.1.2 say BUT BECACUSE THE HEIRARHCY IS BIGGER FOR 3.6.1.1.2 WE INADVERTANTLY CALL THAT
        if actual_score == max_score:
            best_score = actual_score
            best_genotype = genotype
            num_good_nodes = genotype_details[best_genotype]['good_nodes']
        # if actual score is < max score, but is greater than the current best score,
        # then this is our best genotype for the moment
        elif actual_score < max_score and actual_score > best_score:
            best_score = actual_score
            best_genotype = genotype
            num_good_nodes = genotype_details[best_genotype]['good_nodes']

    return best_genotype, num_good_nodes

def extract_lineage_info(lineage_data, genome_name):

    # first, extract phylo spp information - if the spp is not sonnei then don't proceed
    spp_data = lineage_data['species']
    spp_call = list(spp_data.keys())[0]
    # if spp is unknown, then this is not sonnei, exit this function
    if spp_call == "Unknown":
        out_dict = {'genome':[genome_name], 'species':['not S. sonnei'], 'spp_confidence':['NA'], 'name':['NA'],'final_genotype':['NA'],'num_good_nodes':['NA'],'all_genotype_calls':['NA']}
        out_df = pd.DataFrame(out_dict, columns=['genome', 'genotype', 'name', 'num_good_nodes', 'all_genotype_calls'])
        return out_df, spp_call
    else:
        # if it is sonnei, then get the percentage
        spp_percentage = spp_data["Shigella_sonnei"]["percent_coverage"]
        # if the percentage is <90, then exit this function as it's likely not sonnei
        if spp_percentage < 90:
            out_dict = {'genome':[genome_name], 'species':['not S. sonnei'], 'spp_confidence':['NA'], 'name':['NA'],'final_genotype':['NA'],'num_good_nodes':['NA'],'all_genotype_calls':['NA']}
            out_df = pd.DataFrame(out_dict, columns=['genome', 'genotype', 'name', 'num_good_nodes', 'all_genotype_calls'])
            return out_df, "Unknown"

    # if we are sonnei, then continue
    # set up dictionary for final table output
    lineage_out_dict = {'genome': [genome_name]}

    # dictionary of human readable names and their corresponding lineages
    lineage_name_dict = {'lineage1': 'Lineage I', 'lineage1.1': '-', 'lineage1.2': '-', 'lineage1.3': '-', 'lineage1.4': '-', 'lineage1.5': '-', 'lineage1.5.1': '-', 'lineage1.5.2': '-', 'lineage1.5.3': '-', 'lineage1.6': '-', 'lineage1.6.1': '-', 'lineage1.6.2': '-', 'lineage1.6.3': '-', 'lineage1.6.4': '-', 'lineage2': 'Lineage II', 'lineage2.1': '-', 'lineage2.1.1': '-', 'lineage2.1.2': '-', 'lineage2.1.3': '-', 'lineage2.1.4': '-', 'lineage2.1.5': '-', 'lineage2.1.6': '-', 'lineage2.1.7': '-', 'lineage2.1.8': '-', 'lineage2.2': '-', 'lineage2.3': '-', 'lineage2.4': '-', 'lineage2.4.1': '-', 'lineage2.4.2': '-', 'lineage2.4.3': '-', 'lineage2.5': '-', 'lineage2.5.1': '-', 'lineage2.5.2': '-', 'lineage2.6': '-', 'lineage2.6.1': '-', 'lineage2.6.2': '-', 'lineage2.7': '-', 'lineage2.7.1': '-', 'lineage2.7.2': '-', 'lineage2.7.3': '-', 'lineage2.7.4': '-', 'lineage2.8': '-', 'lineage2.8.1': 'Korea II', 'lineage2.8.2': '-', 'lineage2.9': 'Latin America IIa', 'lineage2.9.1': 'Latin America IIa', 'lineage2.9.2': 'Latin America IIa', 'lineage2.10': '-', 'lineage2.10.1': 'Latin America IIa', 'lineage2.10.2': 'Latin America IIa', 'lineage2.10.3': 'Latin America IIa', 'lineage2.10.4': 'Latin America IIa', 'lineage2.10.5': 'Latin America IIa', 'lineage2.10.6': 'Latin America IIa', 'lineage2.10.7': 'Latin America IIa', 'lineage2.10.8': 'Latin America IIa', 'lineage2.10.9': 'Latin America IIa', 'lineage2.11': 'Latin America IIb', 'lineage2.11.1': 'Latin America IIb', 'lineage2.11.2': 'Latin America IIb', 'lineage2.11.3': 'Latin America IIb', 'lineage2.11.4': 'Latin America IIb', 'lineage2.11.5': 'Latin America IIb', 'lineage2.12': 'Latin America IIb', 'lineage2.12.1': 'Latin America IIb', 'lineage2.12.2': 'Latin America IIb', 'lineage2.12.3': 'Latin America IIb', 'lineage2.12.4': 'Latin America IIb', 'lineage3': 'Lineage III', 'lineage3.1': '-', 'lineage3.2': '-', 'lineage3.3': '-', 'lineage3.4': 'Latin America III', 'lineage3.4.1': 'Latin America III', 'lineage3.4.2': 'Latin America III', 'lineage3.4.3': 'Latin America III', 'lineage3.4.4': 'Latin America III', 'lineage3.4.5': 'Latin America III', 'lineage3.4.6': 'Latin America III', 'lineage3.5': '-', 'lineage3.6': 'Central Asia III', 'lineage3.6.1': 'CipR_parent', 'lineage3.6.1.1': 'CipR', 'lineage3.6.1.1.1': 'CipR.SEA', 'lineage3.6.1.1.2': 'CipR.MSM5', 'lineage3.6.1.1.3': 'CipR', 'lineage3.6.1.1.3.1': 'CipR.MSM1', 'lineage3.6.2': 'Central Asia III', 'lineage3.6.3': 'Central Asia III', 'lineage3.6.4': 'Central Asia III', 'lineage3.7': 'Global III', 'lineage3.7.1': 'Global III', 'lineage3.7.3': 'Global III', 'lineage3.7.4': 'Global III', 'lineage3.7.5': 'Global III', 'lineage3.7.6': 'Global III', 'lineage3.7.7': 'Global III', 'lineage3.7.8': 'Global III', 'lineage3.7.9': 'Global III', 'lineage3.7.10': 'Global III', 'lineage3.7.11': 'Global III', 'lineage3.7.12': 'Global III', 'lineage3.7.13': 'Global III', 'lineage3.7.14': 'Global III', 'lineage3.7.15': 'Global III', 'lineage3.7.16': 'Global III', 'lineage3.7.17': 'Global III', 'lineage3.7.18': 'Global III', 'lineage3.7.19': 'Global III', 'lineage3.7.20': 'Global III', 'lineage3.7.21': 'Global III', 'lineage3.7.22': 'Global III', 'lineage3.7.23': 'Global III', 'lineage3.7.24': 'Global III', 'lineage3.7.25': 'MSM4', 'lineage3.7.26': 'Global III', 'lineage3.7.27': 'Global III', 'lineage3.7.28': 'Global III', 'lineage3.7.29': 'Global III VN', 'lineage3.7.29.1': 'Global III VN2', 'lineage3.7.29.1.1': 'Global III VN3', 'lineage3.7.29.1.1.1': 'Global III VN3.KH2', 'lineage3.7.29.1.1.2': 'Global III VN4', 'lineage3.7.29.1.2': 'Global III VN2.MSM2', 'lineage3.7.29.1.2.1': 'Global III VN2.MSM2.Aus', 'lineage3.7.29.1.3': 'Global III VN2.Hue', 'lineage3.7.29.1.4': 'Global III VN2.KH1', 'lineage3.7.29.1.4.1': 'Global III VN2.KH1.Aus', 'lineage3.7.30.1': 'Global III Middle East III', 'lineage3.7.30.2': 'Global III Middle East III', 'lineage3.7.30.3': 'Global III Middle East III', 'lineage3.7.30.4': 'Global III Israel III', 'lineage3.7.30.4.1': 'Global III OJC', 'lineage3.7.31': 'Global III', 'lineage3.7.32': 'Global III', 'lineage3.7.33': 'Global III', 'lineage3.8': '-', 'lineage3.9': '-', 'lineage3.10': '-', 'lineage4': 'Lineage IV', 'lineage5': 'Lineage V', 'lineage5.1.1': '-', 'lineage5.1.2': '-', 'lineage5.1.3': '-', 'lineage5.1.4': '-', 'lineage5.1.5': '-', 'lineage5.1.6': '-', 'lineage5.1': '-'}

    # this try/except statement deals with instances where for some reason there is no lineage output
    # in the json file
    try:
        genotype_calls = lineage_data['lineage']['lineage']
    except KeyError:
        out_dict = {'genome':[genome_name], 'species':['S. sonnei'], 'spp_confidence':[spp_percentage], 'name':['NA'],'final_genotype':['uncalled'],'num_good_nodes':['NA'],'all_genotype_calls':['NA']}
        out_df = pd.DataFrame(out_dict, columns=['genome', 'genotype', 'name', 'num_good_nodes', 'all_genotype_calls'])
        return out_df, spp_call
    # if there are no calls, populate with none
    if len(genotype_calls) == 0:
        lineage_out_dict['final_genotype'] = ['none']
        lineage_out_dict['name']  = ['none']
        lineage_out_dict['num_good_nodes'] = ['NA']
        lineage_out_dict['all_genotype_calls']  = ['NA']
    # if there is just one call, list that - but CHECK that all values in heirarchy are good (>0.5)
    # only report up to level in heirarchy where we have a good call
    if len(genotype_calls) == 1:
        # inspect calls for genotype
        best_genotype, num_good_nodes = inspect_calls(lineage_data['lineage']['calls_summary'])
        # extract genotype from lineage thing, add that to the table
        genotype = best_genotype.split('lineage')[1]
        lineage_out_dict['final_genotype'] = [genotype]
        lineage_out_dict['name'] = [lineage_name_dict[best_genotype]]

        # get the nodes and all calls
        lineage_out_dict['num_good_nodes'] = [num_good_nodes]
        lineage_out_dict['all_genotype_calls'] = genotype_calls
    # if there is more than one call, we want to report the best, but also other calls
    elif len(genotype_calls) > 1:
        # get the info for each call
        call_summary = lineage_data['lineage']['calls_summary']
        # work out which lineage has the best call
        best_genotype, num_good_nodes = inspect_calls(call_summary)
        # now write out info
        lineage_out_dict['final_genotype'] = [best_genotype.split('lineage')[1]]
        lineage_out_dict['name'] = [lineage_name_dict[best_genotype]]
        #lineage_nodes = lineage_data['lineage']['calls_summary'][best_genotype]
        lineage_out_dict['num_good_nodes'] = [num_good_nodes]
        lineage_out_dict['all_genotype_calls'] = [';'.join(call_summary.keys())]
    # add species info
    lineage_out_dict['species']=['S. sonnei']
    lineage_out_dict['spp_confidence']=[spp_percentage]
    lineage_out_df = pd.DataFrame(lineage_out_dict, columns=['genome', 'species', 'spp_confidence', 'final_genotype', 'name', 'num_good_nodes', 'all_genotype_calls'])

    return lineage_out_df, spp_call

def main():

    args = get_arguments()
    # list of tables for each result type
    results_tables = []

    # want a single table of outputs

    # read in json files
    for json_file in args.jsons:
        with open(json_file) as f:
            #with open("ERR017671.mykrobe.json") as f:
            myk_result = json.load(f)
        # get genome name (should first and ONLY key at start of json file)
        ## TODO add if statement in case there is more than one key here (not what we expect)
        genome_name = list(myk_result.keys())[0]
        # extract all the data for this genome
        genome_data = myk_result[genome_name]
        # extract the genotyping information
        lineage_data = genome_data["phylogenetics"]
        lineage_table, spp_call = extract_lineage_info(lineage_data, genome_name)
        #all_lineage_tables.append(lineage_table)
        # extract the qrdr information, only if sonnei
        genome_qrdr_table = extract_qrdr_info(genome_data, genome_name, spp_call)
        #all_qrdr_tables.append(genome_qrdr_table)
        # combine the two tables together
        all_info_table = pd.merge(lineage_table, genome_qrdr_table, on="genome", how="inner")
        results_tables.append(all_info_table)


    # concatenate, re-order columns, and write out
    final_results = pd.concat(results_tables, sort=True)
    final_results.to_csv(args.prefix + "_predictResults.tsv", index=False, sep="\t", columns=["genome", "species", "spp_confidence", "final_genotype", "name", "num_qrdr", "parC_S80I", "gyrA_S83L", "gyrA_S83A", "gyrA_D87G", "gyrA_D87N", "gyrA_D87Y", "num_good_nodes", "all_genotype_calls"])

if __name__ == '__main__':
    main()
