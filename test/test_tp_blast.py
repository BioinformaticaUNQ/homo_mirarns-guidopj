import unittest
import json
from ..src.tp import are_args_passed_valid

class TestTP(unittest.TestCase):

    def test_are_args_passed_valid(self):
        args = {
            "target": "new_target"
        }
        argsJSON = json.dumps(args)
        with self.assertRaises(SystemExit) as cm:
            are_args_passed_valid(argsJSON)
