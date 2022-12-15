import os
import argparse
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast.Applications import NcbimakeblastdbCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
import wget

Entrez.email="guidopujadas@gmail.com"

#MK033612.1
#python3 tp_blast.py -gid MK033612.1 -tg vannamei
#python3 tp_blast.py --fasta-file='../seq_test_adn.fasta' --target='chinensis'
#python3 tp_blast.py --fasta-file='../seq_test_arn.fasta' --target='chinensis'

#---agregar dependencias con versiones
#---parametrizar output del resultado obtenido y por default result.json
#---parametrizar la base de datos contra la que hare blast para ADN

#para ADN corro blast contra la base que me dio el user por default non redundant y eso te devuelve los ADN homologos de
#muchas especies, filtro por el target que quiero y con el gen id de esos genes me fijo en la bd de miRNA cuales mirnas coinciden con ese gene-id
#mature.fasta solamente blast contra mirnas
#documentar repo con casos de uso de como correrlo

parser = argparse.ArgumentParser(prog = 'miRNAs translator',
                    description = 'Translates both gene-ids or fasta file to corresponding miRNAs given a specific target',
                    epilog = 'available params above')
parser.add_argument('-fas', '--fasta-file', help='-> a fasta formatted file with an ADN/ARN sequence')
parser.add_argument('-gid', '--gene-id', help=' -> a gene id if fasta file not present')
parser.add_argument('-plim', '--perc-limit', help=' -> identity percentage limit to filter results. If not given 90%% is used as a default')
parser.add_argument('-ev', '--e-value', help=' -> E value limit to filter results. If not given 0.001 is used as a default')
parser.add_argument('-tg', '--target', help=' -> target used to look at the corresponding miRNA')
parser.add_argument('-out', '--output-file', help=' -> choose an specific location, if not given a result.json file will be created')
parser.add_argument('-db', '--database', help=' -> choose an specific DB when providing ADN fasta file. ' +
	'DBs to choose -> ' +
		'18S_fungal_sequences, 28S_fungal_sequences, Betacoronavirus, ITS_RefSeq_Fungi, ITS_eukaryote_sequences, LSU_eukaryote_rRNA, LSU_prokaryote_rRNA, SSU_eukaryote_rRNA, env_nt, env_nr, 16S_ribosomal_RNA, human_genome, landmark, mito, mouse_genome, nr, nt.(00 - 77), pataa, patnt, pdbaa, pdbnt, ref_euk_rep_genomes, ref_prok_rep_genomes, ref_viroids_rep_genomes, ref_viruses_rep_genomes, refseq_select_rna, refseq_select_prot, refseq_protein, refseq_rna, swissprot, tsa_nr, tsa_nt, taxdb'
)

args = parser.parse_args()

is_miRNAs_exists = os.path.exists('../db/mirnest_targets.gz')
if not is_miRNAs_exists:
	print('downloading miRNAs db')
	url = 'http://rhesus.amu.edu.pl/mirnest/copy/downloads/mirnest_targets.gz'
	filename = wget.download(url)
	os.system('mv mirnest_targets.gz ../db')
	os.system('gunzip --keep ../db/mirnest_targets.gz')
	

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

def getChosenDB(db):
	os.system('../ncbi-blast/bin/update_blastdb.pl --decompress ' + db)
	os.system('tar -xvf *{}*.tar.gz'.format(db))
	print(db + '.ndb')
	return db + '.ndb'
	
def check_gene_id_in_mirnas(gene_ids):
	found = False
	with open('../db/mirnest_targets') as result:
		for gid in gene_ids:
			for line in result:
				miRNA_gene_id=line.split("|")[1]
				if miRNA_gene_id == gid:
					print(line)
					found = True
					break
	if not found:
		print('no matches found')

def run_blast():
	should_check_mirnas = False
	results = []
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
			if args.database != None:
				db = getChosenDB(args.database)
				remote=False
				print('ADN, using {} database'.format(db))
			else:
				db = 'nr'
				remote=True
				print('ADN, using nr database')
			should_check_mirnas = True
		else:
			print('ARN, using miRNAs database')
			db_dir='../db/mature.fasta'
			remote=False
			run_make_db = NcbimakeblastdbCommandline(dbtype="nucl",input_file=db_dir, title="miRNA", out="../db/miRNA")
			db="../db/miRNA"
			make_db_cmd = "../ncbi-blast/bin/" + str(run_make_db)
			os.system(make_db_cmd)

	output_file = args.output_file if args.output_file != None else "result.json"
	percentage_limit = int(args.perc_limit) if args.perc_limit != None else 90 #by default we add entries with identity percetage >= 90
	e_value = float(args.e_value) if args.e_value != None else 0.001 #by default we add entries with e value <= 0.001
	run_blastn = NcbiblastnCommandline(query=file_name, db=db, remote=remote, out=output_file, evalue=e_value, outfmt="10 sseqid pident evalue stitle bitscore")
	blast_cmd = "../ncbi-blast/bin/" + str(run_blastn)
	os.system(blast_cmd)

	with open(output_file) as result:
		found = False
		for line in result:
			entry_value_gene_id=line.split(",")[0]
			entry_value_description=line.split(",")[3]
			if args.target in entry_value_description:
				if should_check_mirnas:
					print("found to compare against miRNAs DB")
					results.append(entry_value_gene_id)
					found = True
				print("gene id: ", entry_value_gene_id)
				print("description: ", entry_value_description)
	
	if should_check_mirnas:
		if not found:
			print('nothing found to compare to')
		else:
			print('looking for miRNAs results')
			check_gene_id_in_mirnas(results)

if __name__ == "__main__":
	are_args_passed_valid(args)
	run_blast()
