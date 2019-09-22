"""Ensures that the correct root path is present."""
import os
import sys

def addPath():
  root_dir = sys.path[0]
  for _ in range(2):
    root_dir, _ = os.path.split(root_dir)
  sys.path.insert(0, root_dir)

addPath()
