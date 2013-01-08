import unittest
import sys
import os
import errno

parent_folder = os.path.join(os.path.dirname(__file__), "..")
sys.path.append(parent_folder)

from talutil import OptionObject
from findTAL_old import RunFindTALOldTask

from testutil import ResultTable, TalentTestCase

try:
    os.makedirs('test_output/findTALOld')
except OSError, e:
    if e.errno != errno.EEXIST:
        raise

class FindTALOldResultTable(ResultTable):
    
    tests_version = 1
    
    def undo_column_2(self, dataset):
        
        del dataset['Unique_RE_sites_in_spacer']

    def undo_column_3(self, dataset):
        
        re_index = dataset.headers.index('Unique_RE_sites_in_spacer')
        
        for i, row in enumerate(dataset):
            
            row = list(row)
            
            re_sites = row[re_index].split(" ")
            re_sites.sort()
            row[re_index] = ' '.join(re_sites)
            
            dataset[i] = row
            
    def undo_column_4(self, dataset):
        
        plus_strand_index = dataset.headers.index('Plus strand sequence')
        
        for i, row in enumerate(dataset):
            
            row = list(row)
            plus_strand = row[plus_strand_index].split(" ")
            
            row[plus_strand_index] = ' '.join([plus_strand[1], plus_strand[2].upper(), plus_strand[3]])
            
            dataset[i] = row

class FindTALOldREResultTable(FindTALOldResultTable):
    tests_version = 2
    
class FindTALOldREKnownResultTable(ResultTable):
    
    tests_version = 1
    
    def undo_column_2(self, dataset):
        
        re_index = dataset.headers.index('Unique_RE_sites_in_spacer')
        
        for i, row in enumerate(dataset):
            
            row = list(row)
            
            re_sites = row[re_index].split(" ")
            re_sites.sort()
            row[re_index] = ' '.join(re_sites)
            
            dataset[i] = row
            
class FindTALOldTestSimpleOutput(TalentTestCase):
    
    test_output_dir = "test_output/findTALOld/"
    expected_output_dir = "known_output/findTALOld/"
    
    default_options = {
        'fasta': 'test_input/findTAL_cases_simple.fasta',
        'outpath': '',
        'logFilepath': '/dev/null',
        'min': None,
        'max': None,
        'arraymin': None,
        'arraymax': None,
        'cupstream': 0,
        't1': True,
        'a2': True,
        'tn': True,
        'gn': True,
        'comp': True,
    }

    def test_simple_a2(self):
        
        result_filename = "simple_a2.txt"
        
        options = self.getDefaultOptions(result_filename)
        options['a2'] = False
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldResultTable, result_filename)
        
    def test_simple_gn(self):
        
        result_filename = "simple_gn.txt"
        
        options = self.getDefaultOptions(result_filename)
        options['gn'] = False
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldResultTable, result_filename)
        
    def test_simple_t1(self):
        
        result_filename = "simple_t1.txt"
        
        options = self.getDefaultOptions(result_filename)
        options['t1'] = False
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldResultTable, result_filename)
        
    def test_simple_tn(self):
        
        result_filename = "simple_tn.txt"
        
        options = self.getDefaultOptions(result_filename)
        options['tn'] = False
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldResultTable, result_filename)
        
    def test_simple_comp(self):
        
        result_filename = "simple_comp.txt"
        
        options = self.getDefaultOptions(result_filename)
        options['comp'] = False
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldResultTable, result_filename)
        
    def test_simple_default(self):
        
        result_filename = "simple_default.txt"
        
        options = self.getDefaultOptions(result_filename)
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldResultTable, result_filename)
        
    def test_simple_tn_gn(self):
        
        result_filename = "simple_tn_gn.txt"
        
        options = self.getDefaultOptions(result_filename)
        options['tn'] = False
        options['gn'] = False
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldResultTable, result_filename)
        
    def test_simple_min5_max15(self):
        
        result_filename = "simple_min5_max15.txt"
        
        options = self.getDefaultOptions(result_filename)
        options['min'] = 5
        options['max'] = 15
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldResultTable, result_filename)
        
    def test_simple_arraymin13_arraymax19_comp(self):
        
        result_filename = "simple_arraymin13_arraymax19_comp.txt"
        
        options = self.getDefaultOptions(result_filename)
        options['arraymin'] = 13
        options['arraymax'] = 19
        options['comp'] = False
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldResultTable, result_filename)
        
    def test_simple_arraymin13_arraymax25(self):
        
        result_filename = "simple_arraymin13_arraymax25.txt"
        
        options = self.getDefaultOptions(result_filename)
        options['arraymin'] = 13
        options['arraymax'] = 25
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldResultTable, result_filename)
        
    def test_simple_arraymin15_arraymax10(self):
        
        result_filename = "simple_arraymin15_arraymax10.txt"
        
        options = self.getDefaultOptions(result_filename)
        options['arraymin'] = 15
        options['arraymax'] = 10
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldResultTable, result_filename)
        
    def test_simple_arraymin18_arraymax20_comp(self):
        
        result_filename = "simple_arraymin18_arraymax20_comp.txt"
        
        options = self.getDefaultOptions(result_filename)
        options['arraymin'] = 18
        options['arraymax'] = 20
        options['comp'] = False
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldResultTable, result_filename)

class FindTALOldTestComplexOutput(TalentTestCase):
    
    test_output_dir = "test_output/findTALOld/"
    expected_output_dir = "known_output/findTALOld/"
    
    default_options = {
        'fasta': 'test_input/findTAL_cases_complex.fasta',
        'outpath': '',
        'logFilepath': '/dev/null',
        'min': None,
        'max': None,
        'arraymin': None,
        'arraymax': None,
        'cupstream': 0,
        't1': True,
        'a2': True,
        'tn': True,
        'gn': True,
        'comp': True,
    }
    
    def test_complex_a2_comp_min12_max13(self):
        
        result_filename = "complex_a2_comp_min12_max13.txt"
        
        options = self.getDefaultOptions(result_filename)
        options['a2'] = False
        options['comp'] = False
        options['min'] = 12
        options['max'] = 13
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldResultTable, result_filename)
        
    def test_complex_comp_a2(self):
        
        result_filename = "complex_comp_a2.txt"
        
        options = self.getDefaultOptions(result_filename)
        options['a2'] = False
        options['comp'] = False
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldResultTable, result_filename)
        
    def test_complex_min20_max25_comp(self):
        
        result_filename = "complex_min20_max25_comp.txt"
        
        options = self.getDefaultOptions(result_filename)
        options['comp'] = False
        options['min'] = 20
        options['max'] = 25
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldResultTable, result_filename)
        
    def test_complex_min8_max10_tn(self):
        
        result_filename = "complex_min8_max10_tn.txt"
        
        options = self.getDefaultOptions(result_filename)
        options['tn'] = False
        options['min'] = 8
        options['max'] = 10
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldResultTable, result_filename)
        
    def test_complex_min8_max10_tn_gn(self):
        
        result_filename = "complex_min8_max10_tn_gn.txt"
        
        options = self.getDefaultOptions(result_filename)
        options['tn'] = False
        options['gn'] = False
        options['min'] = 8
        options['max'] = 10
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldResultTable, result_filename)
        
    def test_complex_t1_a2(self):
        
        result_filename = "complex_t1_a2.txt"
        
        options = self.getDefaultOptions(result_filename)
        options['a2'] = False
        options['t1'] = False
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldResultTable, result_filename)
        
    def test_complex_t1_tn_gn(self):
        
        result_filename = "complex_t1_tn_gn.txt"
        
        options = self.getDefaultOptions(result_filename)
        options['t1'] = False
        options['tn'] = False
        options['gn'] = False
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldResultTable, result_filename)
        
    def test_complex_default(self):
        
        result_filename = "complex_default.txt"
        
        options = self.getDefaultOptions(result_filename)
        options['fasta'] = 'test_input/test_genes.fasta'
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldResultTable, result_filename)
        
    def test_complex_all(self):
        
        result_filename = "complex_all.txt"
        
        options = self.getDefaultOptions(result_filename)
        options['t1'] = False
        options['a2'] = False
        options['tn'] = False
        options['gn'] = False
        options['comp'] = False
        options['min'] = 12
        options['max'] = 13
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldResultTable, result_filename)

#@unittest.skip('')
class FindTALOldTestREOutput(TalentTestCase):
    
    test_output_dir = "test_output/findTALOld/"
    expected_output_dir = "known_output/findTALOld/"
    
    default_options = {
        'fasta': 'test_input/findTAL_cases_complex_RE.fasta',
        'outpath': '',
        'logFilepath': '/dev/null',
        'min': None,
        'max': None,
        'arraymin': None,
        'arraymax': None,
        'cupstream': 0,
        't1': True,
        'a2': True,
        'tn': True,
        'gn': True,
        'comp': True,
    }

    def test_complex_RE(self):
        
        result_filename = "complex_RE.txt"
        
        options = self.getDefaultOptions(result_filename)
        
        RunFindTALOldTask(OptionObject(**options))
        
        self.compareResultFiles(FindTALOldREResultTable, result_filename, expected_table_type=FindTALOldREKnownResultTable)
    
if __name__ == '__main__':
    unittest.main()