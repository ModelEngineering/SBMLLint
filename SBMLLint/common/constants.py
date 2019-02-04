"""Constants for SBMLLint."""
from collections import namedtuple
import os

############### TESTS #####################
TEST_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_DIR = os.path.join(TEST_DIR, "tests")
TEST_FILE = os.path.join(TEST_DIR, "test_file.xml")
NUM_REACTIONS = 111
NUM_PARAMETERS = 27
MAX_REACTANTS = 10
NUM_SPECIES = 32

############## FUNCTIONAL CONSTANTS ##############
# DataFrame Columns
VALUE = "value"

#
MOIETY_SEPARATOR = "_"

# Reaction categories
REACTION_1_1 = "reaction_1_1"
REACTION_n_1 = "reaction_n_1"
REACTION_1_n = "reaction_1_n"
REACTION_n_n = "reaction_n_n"
ReactionCategory = collections.namedtupel('ReactionCategory',
    'category predicate')
REACTION_CATEGORIES = [
    ReactionCategory(category=REACTION_1_1,
        predicate=lambda x,y: (x==1) and (y==1)),
    ReactionCategory(category=REACTION_1_n,
        predicate=lambda x,y: (x==1) and (y>1)),
    ReactionCategory(category=REACTION_n_1,
        predicate=lambda x,y: (x>1) and (y==1)),
    ReactionCategory(category=REACTION_n_n,
        predicate=lambda x,y: (x>1) and (y>1)),
    ]
  
  (1, 1): REACTION_1_1,
