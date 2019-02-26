"""Constants for SBMLLint."""
import collections
import os

############### TESTS #####################
TEST_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_DIR = os.path.join(TEST_DIR, "tests")
TEST_FILE = os.path.join(TEST_DIR, "test_file.xml")
TEST_FILE2 = os.path.join(TEST_DIR, "test_file2.xml")
# test_file3: originally curated_017
TEST_FILE3 = os.path.join(TEST_DIR, "test_file3.xml")
TEST_FILE5 = os.path.join(TEST_DIR, "test_file5.antimony")
TEST_FILE4 = os.path.join(TEST_DIR, "test_file4.xml")
# test_file6: originally curated_050
TEST_FILE6 = os.path.join(TEST_DIR, "test_file6.xml")
# test_file7: originally curated_018
TEST_FILE7 = os.path.join(TEST_DIR, "test_file7.xml")
# test_file8: originally curated_004
TEST_FILE8 = os.path.join(TEST_DIR, "test_file8.xml")
# test_file8: originally curated_008
TEST_FILE9 = os.path.join(TEST_DIR, "test_file9.xml")
NUM_REACTIONS = 111
NUM_PARAMETERS = 27
MAX_REACTANTS = 10
NUM_SPECIES = 32

############## FUNCTIONAL CONSTANTS ##############
# DataFrame Columns
VALUE = "value"
MOIETY = "moiety"

#
MOIETY_SEPARATOR = "_"
MOIETY_DOUBLE_SEPARATOR = MOIETY_SEPARATOR + MOIETY_SEPARATOR
KINETICS_SEPARATOR = ";"
ARC_ARROW = "->"

# Reaction categories
REACTION_1_1 = "reaction_1_1"
REACTION_n_1 = "reaction_n_1"
REACTION_1_n = "reaction_1_n"
REACTION_n_n = "reaction_n_n"
REACTION_BOUNDARY = "reaction_boundary"
ReactionCategory = collections.namedtuple('ReactionCategory',
    'category predicate')

# EmptySet in a reaction (ex. curated model 006)
EMPTYSET = "EmptySet"

# Reaction attributes for each arc in MESGraph
REACTION = "reaction"

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
        predicate=lambda x,y: (x==0) or (y==0)),
    ]


# Directories and files
# Where data files are stored by default
DATA_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.dirname(DATA_DIR)
DATA_DIR = os.path.join(DATA_DIR, "data")
