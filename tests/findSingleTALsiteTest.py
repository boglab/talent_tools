import unittest
import sys
import os
import errno

parent_folder = os.path.join(os.path.dirname(__file__), "..")
sys.path.append(parent_folder)

from talutil import OptionObject
from findSingleTALsite import RunFindSingleTALSiteTask

from testutil import ResultTable, TalentTestCase

try:
	os.makedirs('test_output/findSingleTALsite')
except OSError, e:
	if e.errno != errno.EEXIST:
		raise


class FindSingleTALResultTable(ResultTable):
	
	tests_version = 1
			
	def undo_column_2(self, dataset):
		
		plus_strand_index = dataset.headers.index('Plus strand sequence')
		target_index = dataset.headers.index('Target sequence')
		strand_index = dataset.headers.index('Strand')
		
		for i, row in enumerate(dataset):
			
			row = list(row)
			
			target_strand = row[target_index].split(" ")
			row[target_index] = target_strand[1]
			
			if row[strand_index] == "Plus":
				row[plus_strand_index] = row[target_index]
			else:
				plus_strand = row[plus_strand_index].split(" ")
				row[plus_strand_index] = plus_strand[0]
			
			dataset[i] = row

class FindSingleTALTestSimpleOutput(TalentTestCase):
	
	test_output_dir = "test_output/findSingleTALsite/"
	expected_output_dir = "known_output/findSingleTALsite/"
	
	default_options = {
		'fasta': 'test_input/findSingleTALsite_cases_simple.fasta',
		'outpath': '',
		'logFilepath': '/dev/null',
		'arraymin': None,
		'arraymax': None,
		'cupstream': 0,
		't1': True,
		'a2': True,
		'tn': True,
		'gn': True,
		'comp': True,
		'revcomp': False,
	}
	
	def test_simple_a2(self):
		
		result_filename = "simple_a2.txt"
		
		options = self.getDefaultOptions(result_filename)
		options['a2'] = False
		
		RunFindSingleTALSiteTask(OptionObject(**options))
		
		self.compareResultFiles(FindSingleTALResultTable, result_filename)
		
	def test_simple_t1(self):
		
		result_filename = "simple_t1.txt"
		
		options = self.getDefaultOptions(result_filename)
		options['t1'] = False
		
		RunFindSingleTALSiteTask(OptionObject(**options))
		
		self.compareResultFiles(FindSingleTALResultTable, result_filename)
		
	def test_simple_tn(self):
		
		result_filename = "simple_tn.txt"
		
		options = self.getDefaultOptions(result_filename)
		options['tn'] = False
		
		RunFindSingleTALSiteTask(OptionObject(**options))
		
		self.compareResultFiles(FindSingleTALResultTable, result_filename)
		
	def test_simple_comp(self):
		
		result_filename = "simple_comp.txt"
		
		options = self.getDefaultOptions(result_filename)
		options['comp'] = False
		
		RunFindSingleTALSiteTask(OptionObject(**options))
		
		self.compareResultFiles(FindSingleTALResultTable, result_filename)
		
	def test_simple_default(self):
		
		result_filename = "simple_default.txt"
		
		options = self.getDefaultOptions(result_filename)
		
		RunFindSingleTALSiteTask(OptionObject(**options))
		
		self.compareResultFiles(FindSingleTALResultTable, result_filename)
		
	def test_simple_tn_gn(self):
		
		result_filename = "simple_tn_gn.txt"
		
		options = self.getDefaultOptions(result_filename)
		options['tn'] = False
		options['gn'] = False
		
		RunFindSingleTALSiteTask(OptionObject(**options))
		
		self.compareResultFiles(FindSingleTALResultTable, result_filename)
		
	def test_simple_min13_max19(self):
		
		result_filename = "simple_min13_max19.txt"
		
		options = self.getDefaultOptions(result_filename)
		options['arraymin'] = 13
		options['arraymax'] = 19
		
		RunFindSingleTALSiteTask(OptionObject(**options))
		
		self.compareResultFiles(FindSingleTALResultTable, result_filename)
		
	def test_simple_min13_max25(self):
		
		result_filename = "simple_min13_max25.txt"
		
		options = self.getDefaultOptions(result_filename)
		options['arraymin'] = 13
		options['arraymax'] = 25
		
		RunFindSingleTALSiteTask(OptionObject(**options))
		
		self.compareResultFiles(FindSingleTALResultTable, result_filename)
		
	def test_simple_min18_max20_comp(self):
		
		result_filename = "simple_min18_max20_comp.txt"
		
		options = self.getDefaultOptions(result_filename)
		options['arraymin'] = 18
		options['arraymax'] = 20
		options['comp'] = False
		
		RunFindSingleTALSiteTask(OptionObject(**options))
		
		self.compareResultFiles(FindSingleTALResultTable, result_filename)

class FindSingleTALTestComplexOutput(TalentTestCase):
	
	test_output_dir = "test_output/findSingleTALsite/"
	expected_output_dir = "known_output/findSingleTALsite/"
	
	default_options = {
		'fasta': 'test_input/findSingleTALsite_cases_complex.fasta',
		'outpath': '',
		'logFilepath': '/dev/null',
		'arraymin': None,
		'arraymax': None,
		'cupstream': 0,
		't1': True,
		'a2': True,
		'tn': True,
		'gn': True,
		'comp': True,
		'revcomp': False,
	}
	
	def test_complex_a2_tn(self):
		
		result_filename = "complex_a2_tn.txt"
		
		options = self.getDefaultOptions(result_filename)
		options['a2'] = False
		options['tn'] = False
		
		RunFindSingleTALSiteTask(OptionObject(**options))
		
		self.compareResultFiles(FindSingleTALResultTable, result_filename)
		
	def test_complex_t1_a2(self):
		
		result_filename = "complex_t1_a2.txt"
		
		options = self.getDefaultOptions(result_filename)
		options['t1'] = False
		options['a2'] = False
		
		RunFindSingleTALSiteTask(OptionObject(**options))
		
		self.compareResultFiles(FindSingleTALResultTable, result_filename)
		
	def test_complex_t1_comp(self):
		
		result_filename = "complex_t1_comp.txt"
		
		options = self.getDefaultOptions(result_filename)
		options['t1'] = False
		options['comp'] = False
		
		RunFindSingleTALSiteTask(OptionObject(**options))
		
		self.compareResultFiles(FindSingleTALResultTable, result_filename)
		
	def test_complex_t1_revcomp_tn_gn(self):
		
		result_filename = "complex_t1_revcomp_tn_gn.txt"
		
		options = self.getDefaultOptions(result_filename)
		options['t1'] = False
		options['revcomp'] = True
		options['tn'] = False
		options['gn'] = False

		RunFindSingleTALSiteTask(OptionObject(**options))
		
		self.compareResultFiles(FindSingleTALResultTable, result_filename)
		
	def test_complex_comp_a2_revcomp(self):
		
		result_filename = "complex_comp_a2_revcomp.txt"
		
		options = self.getDefaultOptions(result_filename)
		options['a2'] = False
		options['comp'] = False
		options['revcomp'] = True
		
		RunFindSingleTALSiteTask(OptionObject(**options))
		
		self.compareResultFiles(FindSingleTALResultTable, result_filename)
		
	def test_complex_comp_revcomp(self):
		
		result_filename = "complex_comp_revcomp.txt"
		
		options = self.getDefaultOptions(result_filename)
		options['comp'] = False
		options['revcomp'] = True
		
		RunFindSingleTALSiteTask(OptionObject(**options))
		
		self.compareResultFiles(FindSingleTALResultTable, result_filename)
		
	def test_complex_comp_revcomp_a2_tn_gn(self):
		
		result_filename = "complex_comp_revcomp_a2_tn_gn.txt"
		
		options = self.getDefaultOptions(result_filename)
		options['comp'] = False
		options['revcomp'] = True
		options['a2'] = False
		options['tn'] = False
		options['gn'] = False
		
		RunFindSingleTALSiteTask(OptionObject(**options))
		
		self.compareResultFiles(FindSingleTALResultTable, result_filename)
		
	def test_complex_comp_t1_a2_tn(self):
		
		result_filename = "complex_comp_t1_a2_tn.txt"
		
		options = self.getDefaultOptions(result_filename)
		options['a2'] = False
		options['tn'] = False
		options['t1'] = False
		options['comp'] = False
		
		RunFindSingleTALSiteTask(OptionObject(**options))
		
		self.compareResultFiles(FindSingleTALResultTable, result_filename)
		
	def test_complex_comp_tn_gn(self):
		
		result_filename = "complex_comp_tn_gn.txt"
		
		options = self.getDefaultOptions(result_filename)
		options['gn'] = False
		options['tn'] = False
		options['comp'] = False
		
		RunFindSingleTALSiteTask(OptionObject(**options))
		
		self.compareResultFiles(FindSingleTALResultTable, result_filename)
		

if __name__ == '__main__':
    unittest.main()