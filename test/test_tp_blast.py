import unittest
import json
import sys # added!
sys.path.append("..") # added!
from src.tp_blast import are_args_passed_valid, isARN, isADN
import argparse

class TestTP(unittest.TestCase):

	def test_are_not_args_passed_valid_target(self):
		args = argparse.Namespace()
		args.target = "new_target"
		args.fasta_file = None
		args.gene_id = None

		with self.assertRaises(SystemExit) as cm:
			are_args_passed_valid(args)
	
	def test_are_args_passed_valid_fasta_file(self):
		args = argparse.Namespace()
		args.target = "new_target"
		args.fasta_file = "fasta file"
		args.gene_id = None

		isValid = are_args_passed_valid(args)
		self.assertTrue(isValid)
		
	def test_are_args_passed_valid_fasta_file_gene_id(self):
		args = argparse.Namespace()
		args.target = "new_target"
		args.fasta_file = None
		args.gene_id = "id"

		isValid = are_args_passed_valid(args)
		self.assertTrue(isValid)
	
	def test_isADN(self):

		sequence = "ATGC"
		isValid = isADN(sequence)
		self.assertTrue(isValid)
	
	def test_isARN(self):

		sequence = "AUGC"
		isValid = isARN(sequence)
		self.assertTrue(isValid)
            
if __name__ == '__main__':
    unittest.main()
