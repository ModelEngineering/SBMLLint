from setuptools import setup, find_packages

INSTALL_FULL = [
    "matplotlib",
    "networkx",
    "nose",
    "numpy",
    "pandas",
    "python-libsbml",
    "pyyaml",
    "scipy",
    "tellurium",
    "urllib3",
    "xlrd",
    ]
INSTALL_PARTIAL = list(INSTALL_FULL)
INSTALL_PARTIAL.remove("tellurium")

def sbmllint_setup(install_requires):
  setup(
      name='SBMLLint',
      version='1.0.3',
      author='Woosub Shin, Joseph L. Hellerstein',
      author_email='jlheller@uw.edu',
      packages=find_packages(exclude=['tests', 'analysis',
          'notebook', 'docs']),
      scripts=[
          'SBMLLint/tools/games', 'SBMLLint/tools/moiety_analysis'],
      url='https://github.com/ModelEngineering/SBMLLint',
      description='Linter for SBML models.',
      long_description=open('README.md').read(),
      long_description_content_type='text/markdown',
      package_dir={'SBMLLint': 'SBMLLint'},
      install_requires=install_requires,
      include_package_data=True,
      data_files=[('data/biomodels',
          ['data/biomodels/biomodels.zip'])],
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Developers',      # Define that your audience are developers
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'License :: OSI Approved :: MIT License',   # Again, pick a license
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
        ],
      )

try:
  sbmllint_setup(INSTALL_FULL)
  print("***Full install succeeded. You can use antimony files.")
except Exception as exp:
  print("***Full install failed. Trying partial install without Tellurium.")
  try:
    sbmllint_setup(INSTALL_PARTIAL)
    print("***Partial install succeeded: could not install tellurium.") 
    print("   You cannot use antimony files.")
  except:
    pass
