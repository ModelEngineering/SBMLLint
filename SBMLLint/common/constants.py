"""Constants for Project."""
import collections
import os

PROJECT_NAME = "SBMLLint"

MOIETY_ANALYSIS = "moiety_analysis"
GAMES = "games"

############### COLUMN NAMES ##############
FILENAME = "filename"
IS_STRUCTURED = "is_structured"
NUM_BOUNDARY_REACTIONS = "num_boundary_reactions"
TOTAL_REACTIONS = "total_reactions"
NUM_IMBALANCED_REACTIONS = "num_imbalanced_reactions"
NUM_BALANCED_REACTIONS = "num_balanced_reactions"
FRAC_BALANCED_REACTIONS = "frac_balanced_reactions"
FRAC_BOUNDARY_REACTIONS = "frac_Boundary_reactions"

################ DIRECTORIES #################
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
for _ in range(2):
  PROJECT_DIR = os.path.dirname(PROJECT_DIR)
CODE_DIR = os.path.join(PROJECT_DIR, PROJECT_NAME)
TEST_DIR = os.path.join(PROJECT_DIR, "tests")

################ DATA DIRECTORIES #################
BIOMODELS_DIR = os.path.join(PROJECT_DIR, "data/biomodels")
BIOMODELS_ZIP_FILENAME = "biomodels.zip"
BIGG_DIR = os.path.join(PROJECT_DIR, "data/bigg")
ANALYSIS_MOIETY_ANALYSIS_DIR = os.path.join(PROJECT_DIR, "analysis")
ANALYSIS_MOIETY_ANALYSIS_DIR = os.path.join(
    ANALYSIS_MOIETY_ANALYSIS_DIR, "moiety_analysis")

############### TEST FILES #####################
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
# test_file10: originally "BIOMD0000000231_url.xml"
TEST_FILE10 = os.path.join(TEST_DIR, "test_file10.xml")
# test_file10: originally "BIOMD0000000253_url.xml"
TEST_FILE11 = os.path.join(TEST_DIR, "test_file11.xml")
# test_file12: originally "BIOMD0000000281_url.xml"
TEST_FILE12 = os.path.join(TEST_DIR, "test_file12.xml")
# test_file13: originally "BIOMD0000000035_url.xml"
TEST_FILE13 = os.path.join(TEST_DIR, "test_file13.xml")
# test_file_games_pp1: originally "BIOMD0000000383_url.xml"
TEST_FILE_GAMES_PP1 = os.path.join(TEST_DIR, "test_file_games_pp1.xml")
# test_file_games_pp2: originally "BIOMD0000000018_url.xml"
TEST_FILE_GAMES_PP2 = os.path.join(TEST_DIR, "test_file_games_pp2.xml")
# test_file_games_report1: originally "BIOMD0000000248_url.xml"
TEST_FILE_GAMESREPORT1 = os.path.join(TEST_DIR, "test_file_games_report1.xml")
# test_file_games_report2: originally "BIOMD0000000007_url.xml"
TEST_FILE_GAMESREPORT2 = os.path.join(TEST_DIR, "test_file_games_report2.xml")
# test_file_games_report3: originally "BIOMD0000000167_url.xml"
TEST_FILE_GAMESREPORT3 = os.path.join(TEST_DIR, "test_file_games_report3.xml")
NUM_REACTIONS = 111
NUM_PARAMETERS = 27
MAX_REACTANTS = 10
NUM_SPECIES = 32

############## FUNCTIONAL CONSTANTS ##############
# DataFrame Columns
VALUE = "value"
MOIETY = "moiety"

#
ARC_ARROW = "->"
EQUAL = "="
LABEL_SEPARATOR = ":"
LESSTHAN = "<"
KINETICS_SEPARATOR = ";"
MOIETY_SEPARATOR = "_"
MOIETY_DOUBLE_SEPARATOR = MOIETY_SEPARATOR + MOIETY_SEPARATOR

# Reaction categories
REACTION_1_1 = "reaction_1_1"
REACTION_n_1 = "reaction_n_1"
REACTION_1_n = "reaction_1_n"
REACTION_n_n = "reaction_n_n"
REACTION_BOUNDARY = "reaction_boundary"
# for reduced reactions
REACTION_REDUNDANT = "reaction_redundant"
REACTION_ERROR = "reaction_error"
ReactionCategory = collections.namedtuple('ReactionCategory',
    'category predicate')
ReactionComponents = collections.namedtuple('ReactionComponents',
    'label reactants products')

# EmptySet in a reaction (ex. curated model 006)
EMPTYSET = "EmptySet"

# Reaction attributes for each arc in MESGraph
REACTION = "reaction"

# The selected category is the first one that first has
# a satisfied predicate
# Arguments:
# x: number of reactants, y: number of products
# z: sum of reactants stoichiometry
# w: sum of products stoichiometry
REACTION_CATEGORIES = [
    ReactionCategory(category=REACTION_1_1,
        predicate=lambda x,y,z,w: (x==1) and (y==1) and (z==w)),
    ReactionCategory(category=REACTION_1_n,
        predicate=lambda x,y,z,w: ((x==1) and (y>1) and (z==1.00)) \
                               or ((x==1) and (y==1) and (z<w))),
    ReactionCategory(category=REACTION_n_1,
        predicate=lambda x,y,z,w: ((x>1) and (y==1) and (w==1.00)) \
                               or ((x==1) and (y==1) and (z>w))),
    ReactionCategory(category=REACTION_n_n,
        predicate=lambda x,y,z,w: ((x>1) and (y>1)) 
                               or ((x==1) and (y>1) and (z!=1.00)) \
                               or ((x>1) and (y==1) and (w!=1.00))),
    ReactionCategory(category=REACTION_BOUNDARY,
        predicate=lambda x,y,z,w: (x==0) or (y==0) or (z==0) or (w==0)),
    ]
# Similar to REACTION_CATEGORIES but have different arguments r, p 
# and two additional categories - reaction_redundant and reaction_error
# r: list of reactants stoichiometry
# p: list of products stoichiometry
REACTION_SUMMARY_CATEGORIES = [
    ReactionCategory(category=REACTION_REDUNDANT,
        predicate=lambda x,y,r,p: (x==0) and (y==0)),
    ReactionCategory(category=REACTION_ERROR,
        predicate=lambda x,y,r,p: ((x==0) and (y!=0)) \
                                  or ((x!=0) and (y==0))),
    ReactionCategory(category=REACTION_1_1,
        predicate=lambda x,y,r,p: (x==1) and (y==1) and (sum(r)==sum(p))),
    
    ReactionCategory(category=REACTION_1_n,
        predicate=lambda x,y,r,p: (x==1) and (sum([r[0]<=e for e in p])==len(p))),
                     
    ReactionCategory(category=REACTION_n_1,
        predicate=lambda x,y,r,p: (y==1) and (sum([p[0]<=e for e in r])==len(r))),
    ]


# summarized version of a reaction
ReactionSummary = collections.namedtuple('ReactionSummary', 
                  'label reactants products category')


# building block for each MESGraph path
PathComponents = collections.namedtuple('PathComponents',
                                        'node1 node2 reactions')

# used for creating a report
NULL_STR = ""

# Top level keys in configuration file
CFG_IGNORED_MOIETIES = "ignored_moieties"
CFG_IGNORED_MOLECULES = "ignored_molecules"
CFG_MOIETY_STRUCTURE = "moiety_structure"
CFG_PROCESS_BOUNDARY_REACTIONS = "process_boundary_reactions"
CFG_GAMES_THRESHOLD = "games_threshold_num_reactions"
CFG_SECTIONS = [
    CFG_IGNORED_MOLECULES,
    CFG_IGNORED_MOIETIES,
    CFG_PROCESS_BOUNDARY_REACTIONS,
    CFG_MOIETY_STRUCTURE,
    CFG_GAMES_THRESHOLD,
    ]

# Default values for configuration file
CFG_DEFAULTS = {}
CFG_DEFAULTS[CFG_IGNORED_MOIETIES] = ['DUMMYMOIETY']
CFG_DEFAULTS[CFG_IGNORED_MOLECULES] = ['DUMMYMOLECULE']
CFG_DEFAULTS[CFG_PROCESS_BOUNDARY_REACTIONS] = False
CFG_DEFAULTS[CFG_GAMES_THRESHOLD] = 20
CFG_DEFAULT_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CFG_DEFAULT_PATH = os.path.join(CFG_DEFAULT_PATH, ".sbmllint_cfg.yml")
