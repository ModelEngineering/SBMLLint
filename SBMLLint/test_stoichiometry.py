####
# Test stoichiometry_matrix.py
####


import stoichiometry_matrix as sm

from nose.tools import assert_equals
from nose import with_setup 
import numpy as np
import os
import pandas as pd
import tellurium as te
import tesbml

# To test the file, run the command below on terminal
# nosetests test_stoichiometry.py --with-coverage --cover-package=test_stoichiometry


current_dir = os.getcwd()
data_dir = os.path.abspath(os.path.join(current_dir, os.pardir, 'curated_data'))

# load a model and prepare a StoichiometryMatrix object
num = 15
format_num = format(num, '03d')
document = tesbml.readSBML(os.path.join(data_dir, 'curated_' + format_num + '.xml'))
model = document.getModel()
stoi_class = sm.StoichiometryMatrix(model)
stoi_matrix = stoi_class.buildMatrix()
default_consistency = stoi_class.consistent

# before using isConsistent, default should be None
def test_default_consistency():
	assert_equals(default_consistency, None)

class testSM:

	# type should be pandas DataFrame
	def test_data_type(self):
		assert_equals(type(stoi_matrix), pd.core.frame.DataFrame)

	# number of reaction should match that of model.getNumReactions()
	def test_num_reactions(self):
		assert_equals(len(stoi_matrix.columns), model.getNumReactions())

	# should show consistency test result
	def test_consistency(self):
		assert_equals(stoi_class.isConsistent(), False)


"""
# a model which raises a type error
null_num = 649
null_format_num = format(num, '03d')
null_document = tesbml.readSBML(os.path.join(data_dir, 'curated_' + format_num + '.xml'))
null_model = null_document.getModel()


def test_null():
	assert sm.StoichiometryMatrix(null_model) == TypeError

"""

 
def multiply(x, y):
	return x*y

"""
The following is an example using @with_setup method

def my_setup_function():
    print ("my_setup_function")
 
def my_teardown_function():
    print ("my_teardown_function")

@with_setup(my_setup_function, my_teardown_function)
def test_numbers_3_4():
    print ('test_numbers_3_4  <============================ actual test code')
    assert multiply(3,4) == 12
 
@with_setup(my_setup_function, my_teardown_function)
def test_strings_a_3():
    print ('test_strings_a_3  <============================ actual test code')
    assert multiply('a',3) == 'aaa'


class TestUM:
 
    def setup(self):
        print ("TestUM:setup() before each test method")
 
    def teardown(self):
        print ("TestUM:teardown() after each test method")
 
    @classmethod
    def setup_class(cls):
        print ("setup_class() before any methods in this class")
 
    @classmethod
    def teardown_class(cls):
        print ("teardown_class() after any methods in this class")
 
    def test_numbers_5_6(self):
        print ('test_numbers_5_6()  <============================ actual test code')
        assert multiply(5,6) == 30
 
    def test_strings_b_2(self):
        print ('test_strings_b_2()  <============================ actual test code')
        assert multiply('b',2) == 'bb'
      
 """


