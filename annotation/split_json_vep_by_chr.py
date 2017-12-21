from __future__ import print_function
import argparse
import sys
import json
import os
import pandas as pd
import numpy as np
from collections import OrderedDict
import re
from datetime import datetime
import gzip


def get_csv_row_dict(var_dict, csv_row_keys, csv_json_dict):
    
    csv_row_dict = OrderedDict([(key,np.NaN) for key in csv_row_keys]) 
    csv_row_dict["#Uploaded_variation"] = var_dict["id"].replace("-","_")
    [chr,pos,ref,alt] = var_dict["id"].split("-") 
    [start,end] = sorted([int(x) for x in [var_dict["start"],var_dict["end"]]])
    csv_row_dict["Location"] = "{0}:{1}".format(chr,start) if start == end else "{0}:{1}-{2}".format(chr,start,end)
    allele = var_dict["allele_string"][var_dict["allele_string"].index("/")+1:]
    csv_row_dict["Allele"] = "NA" if allele == "-" else allele
    csv_row_dict["Feature_type"] = "Transcript"
 
    if "transcript_consequences" in var_dict:
        trans_csq_l = var_dict["transcript_consequences"]
        #Sort the transcript consequences.
        trans_csq_l.sort(key=lambda trans_csq: get_transcript_csq_sort_key(trans_csq,var_dict["most_severe_consequence"]))
        csv_row_dict = get_data_from_sub_dict(trans_csq_l[0], csv_row_dict, csv_json_dict)
        if csv_row_dict["CANONICAL"] == 1: csv_row_dict["CANONICAL"] = "YES"
        if "cdna_start" in trans_csq_l[0]:
            csv_row_dict["cDNA_position"] = "{0}".format(trans_csq_l[0]["cdna_start"]) if trans_csq_l[0]["cdna_start"] == trans_csq_l[0]["cdna_end"] else "{0}-{1}".format(trans_csq_l[0]["cdna_start"],trans_csq_l[0]["cdna_end"])

    if "custom_annotations" in var_dict:
        if len(var_dict["custom_annotations"]) > 0:
            if "gnomad_genomes" in var_dict["custom_annotations"]:
                csv_row_dict = get_data_from_sub_dict(var_dict["custom_annotations"]["gnomad_genomes"][0]["fields"], csv_row_dict, csv_json_dict)
            elif "gnomad_exomes" in var_dict["custom_annotations"]:
                csv_row_dict = get_data_from_sub_dict(var_dict["custom_annotations"]["gnomad_exomes"][0]["fields"], csv_row_dict, csv_json_dict)
            if "kaviar" in var_dict["custom_annotations"]:
                csv_row_dict = get_data_from_sub_dict(var_dict["custom_annotations"]["kaviar"][0]["fields"], csv_row_dict, csv_json_dict)
            if "dbsnp" in var_dict["custom_annotations"]:
                csv_row_dict["Existing_variation"] = var_dict["custom_annotations"]["dbsnp"][0]["name"]
                
    return [chr,csv_row_dict]
 
    
def get_transcript_csq_sort_key(transcript_csq,most_severe_csq):
        sort_val = 0
        if most_severe_csq not in transcript_csq["consequence_terms"]:
            sort_val += 1
        if "canonical" in transcript_csq:
            if transcript_csq["canonical"] != 1:
                sort_val += 1
        return sort_val


def get_data_from_sub_dict(sub_dict, csv_row_dict, csv_json_dict):
    for key in csv_row_dict:
        if key in csv_json_dict and csv_json_dict[key] in sub_dict and pd.isnull(csv_row_dict[key]):
            if "[" in str(sub_dict[csv_json_dict[key]]):
	        csv_row_dict[key] = ",".join(sub_dict[csv_json_dict[key]]) 
	    else:
                csv_row_dict[key] = sub_dict[csv_json_dict[key]]        
    return csv_row_dict


def open_new_csv(chr,csv_row_keys):
    log(chr)
    csv = gzip.open(os.path.join(data_dir,"VEP_chr{0}.csv.gz".format(chr)), "w")
    csv.write(",".join(['"{0}"'.format(x) for x in csv_row_keys]) +"\n")
    return csv


def log(txt):
    print("{0}: {1}".format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'),txt)) 


#Set VEP directory.
parser = argparse.ArgumentParser(description="Split VEP json by chr and write to csv.")
parser.add_argument('--uclex_build', type=str, required=True)
args = parser.parse_args() 
data_dir = os.path.join("/SAN/vyplab/UCLex/mainset_{0}".format(args.uclex_build),"mainset_{0}_VEP".format(args.uclex_build))

#Make the columns for the CSV output.
transcript_l = ["#Uploaded_variation","Location","Allele","Consequence","Gene","Feature","Feature_type","cDNA_position","Existing_variation","IMPACT","STRAND","SYMBOL","SYMBOL_SOURCE","HGNC_ID","BIOTYPE","CANONICAL","HGVSc","EXON"]
protein_l = ["CCDS","PROTEIN_ID","SWISSPROT","TREMBL","UNIPARC"]
pathogenicity_l = ["SIFT","PolyPhen","CAROL","CADD_RAW","CADD"]
pop_l = ["AFR","AMR","EAS","FIN","NFE","OTH","SAS"]
exac_l = ["ExAC_MAF"] + ["ExAC_{0}_MAF".format(pop) for pop in pop_l]
gnomad_l = ["gnomad_{0}_{1}".format(stat_type,pop) for stat_type in ["AC","AN","AF","Hom"] for pop in sorted(pop_l + ["ASJ"]) + ["Male","Female","raw"]]
csv_row_keys = transcript_l + protein_l + pathogenicity_l + ["Kaviar"] + exac_l + gnomad_l

#Make a dictionary to map from CSV columns to fields in the json.
transcript_csv_l =  [col for col in transcript_l[transcript_l.index("Gene"):] if col not in ["Feature_type","cDNA_position","Existing_variation"]]
transcript_json_l = [col.lower() for col in transcript_csv_l]
transcript_json_l = [col + "_id" if col=="gene" else "gene_"+col if "symbol" in col else "transcript_id" if col=="feature" else col for col in transcript_json_l]
csv_json_dict = dict(zip(transcript_csv_l,transcript_json_l))
protein_csv_json_dict = dict(zip(protein_l,[col.lower() for col in protein_l]))
path_csv_json_dict = dict(zip(pathogenicity_l,[path_col.lower() + "_prediction" if path_col in ["SIFT","PolyPhen"] else path_col.lower() + "_phred" if path_col=="CADD" else path_col.lower() for path_col in pathogenicity_l]))
exac_csv_json_dict = dict(zip(exac_l,["exac_af" if exac_col.count("_") == 1 else "exac_af_" + exac_col.split("_")[1].lower() for exac_col in exac_l]))
gnomad_csv_json_dict = dict(zip(gnomad_l,[gnomad_col.replace("gnomad_","") for gnomad_col in gnomad_l]))
csv_json_dict.update(protein_csv_json_dict)
csv_json_dict.update(path_csv_json_dict)
csv_json_dict.update(exac_csv_json_dict)
csv_json_dict.update(gnomad_csv_json_dict)
csv_json_dict["Consequence"] = "consequence_terms"
csv_json_dict["Kaviar"] = "AF"

#Read in the json strings and convert to csv output.
current_chr = "1"
csv = open_new_csv(current_chr,csv_row_keys)
csv_data_l = [] 
line_counter = 0
pattern = re.compile("([=:])(NaN)", re.IGNORECASE)
#with open(os.path.join(data_dir,"for_VEP.VEP.head.json")) as infile: #Debugging
with open(os.path.join(data_dir,"for_VEP.VEP.json")) as infile:
    for line in infile:
        line_counter += 1
        line = line.strip()
	d = None
        try:
            line = pattern.sub("\\1null", line)
	    d = json.loads(line)
        except:
            print(line)
        [chr,csv_row_dict] = get_csv_row_dict(d, csv_row_keys, csv_json_dict)
        if chr != current_chr:
            csv.write("\n".join(csv_data_l))
            del csv_data_l[:]
            csv.close()
            current_chr = chr
            csv = open_new_csv(current_chr,csv_row_keys)
        csv_row = ",".join(['"NA"' if pd.isnull(x) else '"{0}"'.format(x) for x in list(csv_row_dict.values())])
	csv_data_l.append(csv_row)
csv.write("\n".join(csv_data_l))
csv.close()

