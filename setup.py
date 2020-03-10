from setuptools import setup, find_packages

INSTALL_REQUIRES = [
    "matplotlib",
    "networkx",
    "nose",
    "numpy",
    "pandas",
    "python-libsbml",
    "pyyaml",
    "scipy",
    "urllib3",
    "xlrd",
    ]

def doSetup(install_requires):
  setup(
      name='SBMLLint',
      version='1.0.11',
      author='Woosub Shin, Joseph L. Hellerstein',
      author_email='jlheller@uw.edu',
      packages=find_packages(exclude=['tests', 'analysis',
          'notebook', 'docs']),
      scripts=[
          'SBMLLint/tools/games',
          'SBMLLint/tools/games.bat',
          'SBMLLint/tools/moiety_analysis',
          'SBMLLint/tools/moiety_analysis.bat',
          'SBMLLint/tools/lp_analysis',
          'SBMLLint/tools/lp_analysis.bat',
          'SBMLLint/tools/make_moiety_structure',
          'SBMLLint/tools/make_moiety_structure.bat',
          'SBMLLint/tools/print_reactions',
          'SBMLLint/tools/print_reactions.bat',
          ],
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


if __name__ == '__main__':
  doSetup(INSTALL_REQUIRES)
