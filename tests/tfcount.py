import sys
import os
import unittest

parent_folder = os.path.join(os.path.dirname(__file__), "..")
sys.path.append(parent_folder)

from talconfig import GENOME_FILE, PROMOTEROME_FILE, BASE_DIR
from btfcount import TargetFinderCountTask

import pickle

#organisms = ["arabidopsis_thaliana", "drosophila_melanogaster", "oryza_sativa"]
organisms = ["oryza_sativa"]

cutoff_values = [3.0, 3.5, 4.0]

rvd_pair_sets = [
    [
        ["HD HD NG HD HD NI HD HD HD NI NI HD NN NN NG NN NG HD NG", "HD HD NI NN NN NI HD NG NN NG NG NG NN HD NG NG"], 
        ["HD HD NG HD HD NI NG NN NG NG NN NI NI HD NI NG", "HD NN NG NG HD NI NI NI NI HD NN HD NG NN NI NG"], 
        ["HD HD NG HD NI NG HD HD NI NN NG NN NI HD NG HD NG", "NI NN HD NI HD NI HD NN HD HD HD NI NG NG NG"], 
    ]
]

min_values = [12, 15, 16]
max_values = [20, 24, 31]

upstream_values = [0, 1]

class CountPairedRvdTALTest(unittest.TestCase):
    pass

def test_generator(test_name, input_filepath, c_upstream, cutoff, rvd_pair_num, rvd_pair, min_value, max_value, known_result):
    
    def test(self):
        
        count_result = TargetFinderCountTask(input_filepath, c_upstream, cutoff, min_value, max_value, [ [ rvd_pair[0].upper(), rvd_pair[1].upper() ] ])
        
        self.assertEqual(count_result[0], known_result)

    return test

if __name__ == '__main__':
       
    for organism in organisms:
    
        with open(BASE_DIR + '/talent/tests/known_output/tfcount/' + organism, 'rb') as input_file:
            offtarget_counts = pickle.load(input_file)
        
        input_filepath = GENOME_FILE % organism
        
        for upstream in upstream_values:
            for i in range(len(rvd_pair_sets)):
                for pair in range(len(rvd_pair_sets[i])):
                    for cutoff_value in cutoff_values:
                        for min_value in min_values:
                            for max_value in max_values:
                                
                                test_key = '_'.join([
                                    'upstream' + str(upstream),
                                    'set' + str(i),
                                    'pair' + str(pair),
                                    'cutoff' + str(cutoff_value).replace('.', '_'),
                                    'min' + str(min_value),
                                    'max' + str(max_value)
                                ])
                                
                                if test_key in offtarget_counts:
                                    
                                    test_name = "%s_%s" % (organism, test_key)
                                    
                                    test = test_generator(test_name, input_filepath, upstream, cutoff_value, i, rvd_pair_sets[i][pair], min_value, max_value, offtarget_counts[test_key])
                                    
                                    setattr(CountPairedRvdTALTest, "test_" + test_name, test)
    
    unittest.main()