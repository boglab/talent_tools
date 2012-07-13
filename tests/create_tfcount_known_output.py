import sys
import os

parent_folder = os.path.join(os.path.dirname(__file__), "..")
sys.path.append(parent_folder)

from talconfig import GENOME_FILE, PROMOTEROME_FILE, BASE_DIR
from paired_talesf import ScorePairedTalesfTask

import pickle

organisms = ["arabidopsis_thaliana", "drosophila_melanogaster", "oryza_sativa"]

cutoff_values = [3.0, 3.5, 4.0]
    [
        ["HD HD NG HD HD NI HD HD HD NI NI HD NN NN NG NN NG HD NG", "HD HD NI NN NN NI HD NG NN NG NG NG NN HD NG NG"], 
        ["HD HD NG HD HD NI NG NN NG NG NN NI NI HD NI NG", "HD NN NG NG HD NI NI NI NI HD NN HD NG NN NI NG"], 
        ["HD HD NG HD NI NG HD HD NI NN NG NN NI HD NG HD NG", "NI NN HD NI HD NI HD NN HD HD HD NI NG NG NG"], 
    ]
]

min_values = [12, 15, 16]
max_values = [20, 24, 31]

upstream_values = [0, 1]

for organism in organisms:
    
    input_filepath = GENOME_FILE % organism
    
    offtarget_counts = {}
    
    for upstream in upstream_values:
        for i in range(len(rvd_pair_sets)):
            for pair in range(len(rvd_pair_sets[i])):
                for cutoff_value in cutoff_values:
                    for min_value in min_values:
                        for max_value in max_values:
                            
                            output_filename = '_'.join([
                                'upstream' + str(upstream),
                                'set' + str(i),
                                'pair' + str(pair),
                                'cutoff' + str(cutoff_value).replace('.', '_'),
                                'min' + str(min_value),
                                'max' + str(max_value)
                            ])
                            
                            output_filepath = BASE_DIR + '/talent/tests/known_output/tfcount/full/' + output_filename
                            
                            ScorePairedTalesfTask(input_filepath, rvd_pair_sets[i][pair][0].upper(), rvd_pair_sets[i][pair][1].upper(), output_filepath, "NA", upstream, cutoff_value, min_value, max_value, 4, organism)
                            
                            if os.path.exists(output_filepath + '.txt'):
                                
                                with open(output_filepath + '.txt', 'rb') as result_file:
                                    
                                    result_file.readline()
                                    result_file.readline()
                                    result_file.readline()
                                    
                                    raw_row = result_file.readline()
                                    
                                    totals = [0, 0, 0, 0, 0]
                                
                                    while raw_row:
                                        
                                        row = raw_row.split('\t')
                                        
                                        first = (1 if row[1] == 'RVD2' else 0)
                                        second = (1 if row[2] == 'RVD2' else 0)
                                        
                                        index = (2 * first + second)
                                        
                                        totals[index] += 1
                                        
                                        raw_row = result_file.readline()
                                    
                                    for j in range(4):
                                        totals[4] += totals[j]
                                    
                                    print(totals)
                                    
                                    offtarget_counts[output_filename] = totals
                                
                                os.remove(output_filepath + '.txt')
                                
                                try:
                                    os.remove(output_filepath + '.gff')
                                except OSError:
                                    try:
                                        os.remove(output_filepath + '.bed')
                                    except OSError:
                                        pass
                                
                                os.remove(output_filepath + '.gff3')
    
    with open(BASE_DIR + '/talent/tests/known_output/tfcount/' + organism, 'wb') as output_file:
        pickle.dump(offtarget_counts, output_file)
