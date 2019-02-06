"""Constants for SBMLLint."""
import collections
import os

############### TESTS #####################
TEST_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_DIR = os.path.join(TEST_DIR, "tests")
TEST_FILE = os.path.join(TEST_DIR, "test_file.xml")
TEST_FILE2 = os.path.join(TEST_DIR, "test2_file.xml")
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
REACTION_BOUNDARY = "reaction_boundary"
ReactionCategory = collections.namedtuple('ReactionCategory',
    'category predicate')

# Molecule-Stoichiometry tuple in Reaction 
MoleculeStoichiometry = collections.namedtuple('MoleculeStoichiometry',
    'molecule stoichiometry')

# EmptySet in a reaction (ex. curated model 006)
EMPTYSET = "EmptySet"

# The selected category is the first one that first has
# a satisfied predicate
REACTION_CATEGORIES = [
    ReactionCategory(category=REACTION_1_1,
        predicate=lambda x,y: (x==1) and (y==1)),
    ReactionCategory(category=REACTION_1_n,
        predicate=lambda x,y: (x==1) and (y>1)),
    ReactionCategory(category=REACTION_n_1,
        predicate=lambda x,y: (x>1) and (y==1)),
    ReactionCategory(category=REACTION_n_n,
        predicate=lambda x,y: (x>1) and (y>1)),
    ReactionCategory(category=REACTION_BOUNDARY,
        predicate=lambda x,y: (x==0) or (y==0))
    ]


# Directories and files
# Where data files are stored by default
DATA_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.dirname(DATA_DIR)
DATA_DIR = os.path.join(DATA_DIR, "data")
