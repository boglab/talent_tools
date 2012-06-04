import tablib
import unittest

class ResultTable(object):
	
	tests_version = 1
	
	output_sorted = False
	
	def get_undo_columns(self):
		return sorted([int(m[12:]) for m in dir(self) if m.startswith("undo_column_")])
	
	def get_undo_comments(self):
		return sorted([int(m[13:]) for m in dir(self) if m.startswith("undo_comment_")])

	def get_comparable(self, results_filepath):
		
		rows = []
		
		with open(results_filepath, "r") as results_file:
			
			results_line = results_file.readline()
			
			while results_line:
				rows.append(results_line.rstrip())
				results_line = results_file.readline()
		
		# This would be the place to apply any undo_ that needs the unsplit rows
		
		while len(rows) > 0 and len(rows[0].split("\t")) == 1:
			rows.pop(0)
		
		headers = rows.pop(0).split("\t")
		
		if not self.output_sorted:
			rows.sort()
		
		split_rows = []
		
		for row in rows:
			split_rows.append(row.split("\t"))
			
		dataset = tablib.Dataset(*split_rows, headers=headers)
		
		# Apply column revisions
		column_revisions = self.get_undo_columns()
		
		while True:
			
			if len(column_revisions) > 0:
				revision = column_revisions.pop()
			else:
				break
			
			if revision > self.tests_version:
				revision_function = getattr(self, "undo_column_" + str(revision))
				revision_function(dataset)
			else:
				break
			
		return dataset

class SortedResultTable(ResultTable):
	output_sorted = True
	
class TalentTestCase(unittest.TestCase):
	
	def assertResultTablesEqual(self, first_table, second_table):
		
		for i in range(len(first_table)):
			self.assertEqual(first_table[i], second_table[i])
	
	def compareResultFiles(self, table_type, result_filename, expected_table_type=ResultTable):
		
		test_result = table_type()
		test_result_table = test_result.get_comparable(self.test_output_dir + result_filename)
		
		expected_result = expected_table_type()
		expected_result_table = expected_result.get_comparable(self.expected_output_dir + result_filename)
		
		self.assertResultTablesEqual(expected_result_table, test_result_table)
		
	def getDefaultOptions(self, result_filename):
		
		options = self.default_options.copy()
		if 'outpath' in options:
			options['outpath'] = (self.test_output_dir + result_filename)
		else:
			options['outputFilepath'] = (self.test_output_dir + result_filename)
			
		return options


