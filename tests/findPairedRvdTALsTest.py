import unittest
import sys
import os
import errno

parent_folder = os.path.join(os.path.dirname(__file__), "..")
sys.path.append(parent_folder)

from talutil import OptionObject
from findPairedRvdTALs import RunPairedTalesfTask

from testutil import ResultTable, SortedResultTable, TalentTestCase

try:
    os.makedirs('test_output/findPairedRvdTALs')
except OSError, e:
    if e.errno != errno.EEXIST:
        raise

class FindPairedRvdTALResultTable(ResultTable):
    
    tests_version = 1
    
    output_sorted = True

    def undo_column_2(self, dataset):
        
        tal1_target_index = dataset.headers.index('TAL 1 Target')
        tal2_target_index = dataset.headers.index('TAL 2 Target')
        
        for i, row in enumerate(dataset):
            
            row = list(row)
            
            tal1_target_strand = row[tal1_target_index].split(" ")
            row[tal1_target_index] = tal1_target_strand[1] if len(tal1_target_strand[0]) == 1 else tal1_target_strand[0]
            
            tal2_target_strand = row[tal2_target_index].split(" ")
            row[tal2_target_index] = tal2_target_strand[1] if len(tal2_target_strand[0]) == 1 else tal2_target_strand[0]
            
            dataset[i] = row

class FindPairedRvdTALTestOutput(TalentTestCase):
    
    test_output_dir = "test_output/findPairedRvdTALs/"
    expected_output_dir = "known_output/findPairedRvdTALs/"
    
    default_options = {
        'fasta': 'test_input/findPairedRvdTALs.fasta',
        'outputFilepath': '',
        'logFilepath': '/dev/null',
        'min': None,
        'max': None,
        'cupstream': 0,
        'genome': False,
        'promoterome': False,
        'cutoff': 3.0,
        'rvdString': '',
        'rvdString2': '',
        'organism': '',
    }

def test_generator(test_name, cutoff, rvd_pair_num, rvd_pair, min_value, max_value):
    
    def test(self):
        
        options = self.getDefaultOptions(test_name)
        
        options['cutoff'] = cutoff
        options['rvdString'] = rvd_pair[0]
        options['rvdString2'] = rvd_pair[1]
        options['min'] = min_value
        options['max'] = max_value
        
        RunPairedTalesfTask(OptionObject(**options))
        
        result_filename = test_name + ".txt"
        
        self.compareResultFiles(FindPairedRvdTALResultTable, result_filename, expected_table_type = SortedResultTable)
        
    return test

if __name__ == '__main__':
    
    cutoff_values = [3.0, 3.5, 4.0]

    rvd_pairs = [
        ["NN HD HD HD NI NG NN HD NG NN NI NN NI NN NN", "NN HD NN NG NN NG NN HD NN NG NN NG NN NG NN NG"],
    ]
    
    min_values = [12, 15, 16]
    max_values = [20, 24, 31]
    
    for cutoff_value in cutoff_values:
        for i in range(len(rvd_pairs)):
            for min_value in min_values:
                for max_value in max_values:
    
                    test_name = '_'.join([
                        'cutoff' + str(cutoff_value).replace('.', '_'),
                        'pair' + str(i),
                        'min' + str(min_value),
                        'max' + str(max_value)
                    ])
                    
                    test = test_generator(test_name, cutoff_value, i, rvd_pairs[i], min_value, max_value)
                    
                    setattr(FindPairedRvdTALTestOutput, "test_" + test_name, test)
    
    unittest.main()
