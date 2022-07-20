#!/usr/bin/env python# 
#-*- coding: utf-8 -*

import sys
import pandas as pd
from Bio import SeqIO
import pyodbc
import paramiko
import time
from datetime import datetime
from regions_dict import region_dict, at_states


'''
script takes bulk downloaded data from GISAID (sequences.fasta and metadata.tsv), filters them according to selected criteria
and outputs filtered sequences to standard output and nextstrain-formatted metadata to generated_metadata.tsv

usage: python3 filter_GISAID.py > output.fasta
'''

"""
Inputs and outputs
"""
# I/O names
input_fasta = "sequences_clean.fasta.bgz"
input_metadata = "metadata.tsv"
output_unformatted_metadata = "root_extended.tsv"
output_metadata = "root_formatted_extended.tsv"

def read_local_sequences(db_accessions):
	# read sequences local
	records_local = SeqIO.index_db(local_fasta+".index", filenames=[local_fasta], format = "fasta")
	local_fasta = "all.fasta.bgz"
	added_seq_local = []
	for ID in db_accesions:
		try:
			if records_local[ID].id not in added_seq_local:
				print(">"+records_local[ID].id)
				print(records_local[ID].seq)
				added_seq_local.append(records_local[ID].id)
			else:
				sys.stderr.write("Ignoring duplicated record: "+ID+"\n")
		except Exception as e:
			sys.stderr.write(e)
			sys.stderr.write("ERROR: Cannot add sequence: "+ID+"\n")
	return added_seq_local


def print_dbase_fastas(dbase_fastas, IDs):
	for ID, fasta in zip(IDs, dbase_fastas):
		print( ">"+ID)
		print (fasta)
	return None


# get database records in pandas df
"""
configure by entering credentials
"""
def get_dbase():
	server = '' 
	database = '' 
	username = '' 
	password = '' 
	conn = pyodbc.connect('DRIVER={ODBC Driver 17 for SQL Server};SERVER='+server+';DATABASE='+database+';UID='+username+';PWD='+ password)
	cursor = conn.cursor()
	df = pd.read_sql_query("SELECT * FROM samples LEFT JOIN seqData ON samples.RNA_ID_int = seqData.RNA_ID_int LEFT JOIN fasta ON seqData.ID = fasta.Sample_ID",conn)
	df = df.loc[:,~df.columns.duplicated()]
	return df


# formatting of the location fields
def break_location(loc_string):
	try:
		loc_string = loc_string.replace(" ","")
		tmp = loc_string.split("/")
		try:
			if tmp[1] in region_dict:
				region = region_dict[tmp[1]]
			else:
				region = tmp[0]
		except:
			region = "N/A"
		try:
			country = tmp[1]
		except:
			country = "N/A"
		try:
			federal_state = tmp[2]
		except:
			federal_state = "N/A"
		try:
			if federal_state == "Vienna":
				location = "Wien-Stadt"
			if country =="Austria" and federal_state not in at_states:
				location = "N/A"
			else:	
				location = tmp[3]
		except:
			location = "N/A"
		
		return region, country, federal_state, location
	except:
		return "N/A","N/A","N/A","N/A"


# Filtering criteria
"""Define your filtering logic here. Logic operators for Pandas are as follows: OR = | , AND = &"""
def filter_dbase(df, pango_lineage= None, country = "Austria", host = "human", AA_changes = None):
	df = df[df["host"]== host]
	df = df[df["sample_date"].isin([None, "N/A","?","",[]])==False]
	df = df[df["status"]=="passed_qc"]
	df = df[df["N_in_Consensus"]<500]
	if pango_lineage != None:
		df = df[df["pangolin_lineage_current"].isin(pango_lineage)]
	if AA_changes != None:
		df = df[df["aminoacid_changes"].str.contains(AA_changes, case = False)]
	df = df[df["sample_source_location_country"] == country]
	return df

def format_aa(string):
	lst = string.replace("(","").replace(")","").split(",")
	return lst

def convert_df(df, records):
	"""format datasets to fit Nextstrain format"""
	def is_avail(id):
		if id in records:
			return True
		else:
			print(id)
			return False
	
	# create sequence ID
	df["Virus name"] =  df["Virus name"].str.replace(' ', '_')
	df.drop_duplicates(subset=["Virus name"], keep=False, inplace = True)
	
	# format metadata data frame according to Nextstrain
	df.loc[:,"strain"] = df["Virus name"]
	df.loc[:,"virus"] = df["Type"]
	df.loc[:,"gisaid_epi_isl"] = df["Accession ID"]
	df.loc[:,"genbank_accession"] = "?"
	df.loc[:,"date"] = df["Collection date"]
	df.loc[:,"region"], df.loc[:,"country"], df.loc[:,"division"], df.loc[:,"location"] = zip(*df["Location"].apply(break_location))
	df.loc[:,"region_exposure"] = "NA"
	df.loc[:,"country_exposure"] = "NA"
	df.loc[:,"division_exposure"] = "NA"
	df.loc[:,"location_exposure"] = "NA"
	df.loc[:,"segment"] = "genome"
	df.loc[:,"length"] = df["Sequence length"]
	df.loc[:,"host"] = df["Host"]
	df.loc[:,"age"] = df["Patient age"]
	df.loc[:,"sex"] = df["Gender"]
	df.loc[:,"Nextstrain_clade"] = "?"
	df.loc[:,"pangolin_lineage"] = df["Pango lineage"]
	df.loc[:,"GISAID_clade"] = "?"
	df.loc[:,"originating_lab"] = "NA"
	df.loc[:,"submitting_lab"] = "NA"
	df.loc[:,"authors"] = "NA"
	df.loc[:,"url"] = "NA"
	df.loc[:,"title"] = "NA"
	df.loc[:,"paper_url"] = "NA"
	df.loc[:,"date_submitted"] = "NA"
	df.loc[:,"purpose_of_sequencing"] = "NA"
	df.to_csv(output_unformatted_metadata, sep = "\t", index = False)
	# check whether sequence is actually present in the fasta file
	df["available"] = df["Virus name"].apply(is_avail)
	df = df[df["available"] == True]
	return df

def convert_dbase_df(df, records = None):
	"""format database records to fit reference dataset"""
	def is_avail(id):
		if id in records:
			return True
		else:
			print(id)
			return False
	# create sequence ID
	df.loc[:,"unified"] =  df["BSF_sample_name"]
	# remove white spaces from unified ID
	df["unified"] = df["unified"].replace(' ', '_', regex=True)
	# remove all duplicates
	df.drop_duplicates(subset=["unified"], keep=False, inplace = True)
	# format metadata data frame according to Nextstrain
	df.loc[:,"strain"] = df["unified"]
	df.loc[:,"virus"] = "betacoronavirus"
	df.loc[:,"gisaid_epi_isl"] = "N/A"
	df.loc[:,"genbank_accession"] = "?"
	df.loc[:,"date"] = df["sample_date"]
	df.loc[:,"region"] = "Austria"
	df.loc[:,"country"] = df["sample_source_location_country"].str.replace(" ","")
	df.loc[:,"division"] = df["sample_source_location_state"].str.replace(" ","")
	df.loc[:,"location"] = df["sample_source_location_district"].str.replace(" ","")
	df.loc[:,"region_exposure"] = "NA"
	df.loc[:,"country_exposure"] = "NA"
	df.loc[:,"division_exposure"] = "NA"
	df.loc[:,"location_exposure"] = "NA"
	df.loc[:,"segment"] = "genome"
	df.loc[:,"length"] = "30000"
	df.loc[:,"host"] = df["host"]
	df.loc[:,"age"] = df["age"]
	df.loc[:,"sex"] = df["sex"]
	df.loc[:,"Nextstrain_clade"] = "?"
	df.loc[:,"pangolin_lineage"] = df["pangolin_lineage_current"]
	df.loc[:,"GISAID_clade"] = "?"
	df.loc[:,"originating_lab"] = df["sample_source_institute"]
	df.loc[:,"submitting_lab"] = "CeMM"
	df.loc[:,"authors"] = "NA"
	df.loc[:,"url"] = "NA"
	df.loc[:,"title"] = "NA"
	df.loc[:,"paper_url"] = "NA"
	df.loc[:,"date_submitted"] = "NA"
	df.loc[:,"purpose_of_sequencing"] = "NA"
	df_sel = df[["strain", "virus", "gisaid_epi_isl", "genbank_accession", "date", "region", "country", "division", "location", "region_exposure", "country_exposure", "division_exposure", "segment", "length", "host", "age", "sex", "Nextstrain_clade", "pangolin_lineage", "GISAID_clade", "originating_lab", "submitting_lab", "authors", "url", "title", "paper_url", "date_submitted", "purpose_of_sequencing"]]
	return df_sel

def subsample(df, n):
	"""Select subsampling fraction for each country in frac= """
	df_master = pd.DataFrame()
	for country in list(df.country.unique()):
		df_sample = df[df["country"]==country]
		if len(df_sample)>n:
			df_sample = df_sample.sample(n=n)
		df_master = pd.concat([df_master, df_sample])
	return df_master
	"""End of subsampling"""
	
def subsample_variants(df, n):
	"""If you want to subsample by country, uncoment following block. Select subsampling fraction for each country in frac= """
	df_master = pd.DataFrame()
	
	for variant in list(df["pangolin_lineage"].unique()):
		df_sample = df[df["pangolin_lineage"]==variant]
		if (len(df_sample)<20) and (variant not in ["B.1.1.7", "P.1", "B.1.617.2", "B.1.351"]):
			continue
		if len(df_sample)>n:
			df_sample = df_sample.sample(n=n)
		
		df_master = pd.concat([df_master, df_sample])
	
	return df_master
	"""End of subsampling"""

def main():
	# read GISAID bulk metadata file
	metadata = pd.read_csv(input_metadata, sep = "\t")
	# read sequences file
	records = SeqIO.index_db(input_fasta+".index", filenames=[input_fasta], format = "fasta")
	
	# build dataset with Austrian sequences
	df_at = metadata[
			
		
				(
					(metadata["Host"].str.contains("human", case = False)) &
					(metadata["Pango lineage"] == "B.1.1.7") &
					(metadata["Is complete?"] == True) &
					(metadata["Is high coverage?"] == True) &
					(metadata["Collection date"].isin([None, "N/A","?","",[]])==False) &
					(metadata["Virus name"].str.contains("CeMM", case = True)==True) &
					(metadata["Collection date"]>"2020-02-01") &
					(metadata["Location"].str.contains("Austria", case = False)==True) #&
					(metadata["AA Substitutions"].str.contains("E484K", case = False)==True) 
			
				) |
				(
					(metadata["Host"].str.contains("human", case = False)) &
					(metadata["Location"].str.contains("Austria", case = False)) &
					(metadata["Is complete?"] == True) &
					(metadata["Is high coverage?"] == True) &
					(metadata["Pango lineage"] == "B.1.1.7") 
			
				 ) 
			]

	# build complementary dataset with European sequences
	df_eu = 	metadata[	
					(metadata["Host"].str.contains("human", case = False)) &
					(metadata["Pango lineage"]=="B.1.1.7") &
					(metadata["Collection date"].isin([None, "N/A","?","",[]])==False) &
					(metadata["Is complete?"] == True) &
					(metadata["Is high coverage?"] == True) &
					(metadata["Virus name"].str.contains("CeMM", case = False)==False) &
					(metadata["Location"].str.contains("Europe", case = False)==True) #&
				] 
	
	# format AT dataset
	df_at = convert_df(df_at, records)
	df_at = df_at.sample(n = 1000)
	# format EU dataset
	df_eu = convert_df(df_eu, records)
	df_eu = df_eu.sample(n= 500)
	
	include = df_at["strain"]
	include.to_csv("include.txt", index=False, header = False)
	df = pd.concat([df_eu,df_aux, df_at])
	
	# remove duplicated records
	df.drop_duplicates(subset=["strain"], keep = "first", inplace = True)
	
	# get list of sequence IDs to select
	sel_accesions = list(df["strain"])
	
	# output metadata dataframe with all fields
	df.to_csv(output_unformatted_metadata, sep = "\t", index = False)
	
	# check if all IDs are unique
	if not df["strain"].is_unique:
		sys.stderr.write("ERROR: duplicated strain IDs: ")
		sys.stderr.write(" ".join(list(df[df.duplicated(["Virus name"])]["Virus name"])))
		quit()
	
	# filter out redundant columns
	df_sel = df[["strain", "virus", "gisaid_epi_isl", "genbank_accession", "date", "region", "country", "division", "location", "region_exposure", "country_exposure", "division_exposure", "segment", "length", "host", "age", "sex", "Nextstrain_clade", "pangolin_lineage", "GISAID_clade", "originating_lab", "submitting_lab", "authors", "url", "title", "paper_url", "date_submitted", "purpose_of_sequencing"]]
	df_sel.to_csv(output_metadata, sep = "\t", index = False)
	
	'''
	Filtering of the sequences file. Will print only sequences selected from metadata.
	'''
	added_seq = []
	
	# Creates an index of the FASTA on the first run.
	#records = SeqIO.index_db(input_fasta+".index", filenames=[input_fasta], format = "fasta")
	for ID in sel_accesions:
		try:
			if records[ID].id not in added_seq:
				print(">"+records[ID].id)
				print(records[ID].seq)
				added_seq.append(records[ID].id)
			else:
				sys.stderr.write("Ignoring duplicated record: "+ID+"\n")
		except Exception as e:
			sys.stderr.write("ERROR: Cannot add sequence: "+ID+"\n")


if __name__ == "__main__":
	main()

