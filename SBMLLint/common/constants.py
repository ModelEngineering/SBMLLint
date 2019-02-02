"""Constants for SBMLLint."""
import os

############### TESTS #####################

TEST_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_DIR = os.path.join(TEST_DIR, "tests")
TEST_FILE = os.path.join(TEST_DIR, "chemotaxis.xml")
NUM_REACTIONS = 111
NUM_PARAMETERS = 27
MAX_REACTANTS = 10
NUM_SPECIES = 32
