from setuptools import setup, find_packages

setup(
    name='SBMLLint',
    version='1.0.0',
    author='Woosub Shin, Joseph L. Hellerstein',
    author_email='jlheller@uw.edu',
    packages=find_packages(),
    scripts=[],
    url='https://github.com/ModelEngineering/SBMLLint',
    description='Linter for SBML models.',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy",
        "pandas",
        "urllib3",
        "networkx",
        "xlrd",
        "nose",
        "pyyaml",
        "libsbml",
        "tellurium",
        ],
    )
