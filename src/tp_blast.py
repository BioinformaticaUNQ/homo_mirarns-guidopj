import os
import argparse
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastnCommandline

Entrez.email="guidopujadas@gmail.com"

#MK033612.1
#python3 tp_blast.py -gid MK033612.1 -tg vannamei
#python3 tp_blast.py --fasta-file='../seq_test_adn.fasta' --target='chinensis'
#python3 tp_blast.py --fasta-file='../seq_test_arn.fasta' --target='chinensis'

parser = argparse.ArgumentParser()
parser.add_argument('-fas', '--fasta-file')
parser.add_argument('-gid', '--gene-id')
parser.add_argument('-plim', '--perc-limit')
parser.add_argument('-ev', '--e-value')
parser.add_argument('-tg', '--target')

args = parser.parse_args()

def isADN(seq):
	return 'T' in seq or 't' in seq
	
def isARN(seq):
	return 'U' in seq or 'u' in seq
	
def get_sequence_by_gene_id(record, gene_id):
	path = ".cache"
	isExist = os.path.exists(path)
	if not isExist:
		os.makedirs(path)
		print(f"The new directory {path} is created!")
	file_name = f".cache/{gene_id}.fasta"
	SeqIO.write(record, file_name, "fasta")
	return file_name

def are_args_passed_valid(args):
	if args.target == None:
		print("You need to provide a system target")
		print("EXIT")
		exit()
	if args.fasta_file == None and args.gene_id == None:
		print("You need to provide either a fasta file or a gene id")
		print("EXIT")
		exit()
	return True

def run_blast():
	if args.gene_id != None:
		handle = Entrez.efetch(db="nucleotide", id=args.gene_id, rettype="gb", retmode="text")
		record = SeqIO.read(handle, "genbank")
		file_name = get_sequence_by_gene_id(record, args.gene_id)
		db='nt'
		remote=True
	else:
		file_name = args.fasta_file
		fasta_parsed_file = list(SeqIO.parse(open(args.fasta_file),'fasta'))
		sequence = fasta_parsed_file[0].seq
		if isADN(sequence):
			print('ADN, using nr database')
			db='nr'
			remote=True
		else:
			print('ARN, using miRNAs database')
			db_dir='../db/mature.fasta'
			remote=False
			run_make_db = NcbimakeblastdbCommandline(dbtype="nucl",input_file=db_dir, title="miRNA", out="../db/miRNA")
			db="miRNA"
			make_db_cmd = "../ncbi-blast/bin/" + str(run_make_db)
			os.system(make_db_cmd)

	percentage_limit = int(args.perc_limit) if args.perc_limit != None else 90 #by default we add entries with identity percetage >= 90
	e_value = int(args.e_value) if args.e_value != None else 0.001 #by default we add entries with e value <= 0.001
	run_blastn = NcbiblastnCommandline(query=file_name, db=db, remote=remote, out="result.json", evalue=e_value, outfmt="10 sseqid pident evalue stitle bitscore")
	blast_cmd = "../ncbi-blast/bin/" + str(run_blastn)
	os.system(blast_cmd)

	with open('result.json') as result:
		for line in result:
			entry_value_description=line.split(",")[3]
			if args.target in entry_value_description:
				print(entry_value_description)

if __name__ == "__main__":
	are_args_passed_valid(args)
	run_blast()
