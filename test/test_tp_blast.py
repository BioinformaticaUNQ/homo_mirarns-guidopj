import unittest
import json
import sys # added!
sys.path.append("..") # added!
from src.tp_blast import are_args_passed_valid

class TestTP(unittest.TestCase):

	def test_are_args_passed_valid(self):
		args = {
			"target": "new_target"
		}

		with self.assertRaises(SystemExit) as cm:
			are_args_passed_valid(args)
            
if __name__ == '__main__':
    unittest.main()
