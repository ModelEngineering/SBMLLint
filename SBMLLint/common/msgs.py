"""Generate messages for the application."""

import sys

def error(text):
  print("***Error in SBMLLint. Reason follows.")
  print("   %s" % text)
  sys.exit()
