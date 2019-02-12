"""Initializes for using SBMLLint."""

import os
import sys


PROJECT_DIR = os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_DIR)
