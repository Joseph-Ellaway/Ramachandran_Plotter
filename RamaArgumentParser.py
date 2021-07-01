<<<<<<< HEAD
''' ================================================================================================================

	If running script from the command line, functions here are called to parsed user's arguments into the main() function
	in RamachandranPlotter.py.
	
	Version 2.0.1:
	 - Relies on the easily accessible Biopython package, rather than Phenix as in versions <2.0
	 - User arguments can be now easily parsed in from the command line (as facilitated by functions here)
	 - If required, the script could be implemented into existing protein analysis pipelines by importing this 
	 function ( main() )

	Author information:
	 - Joseph I. J. Ellaway
	 - josephellaway@gmail.com
	 - https://github.com/Joseph-Ellaway

	================================================================================================================ '''

import argparse



def CollctUserArgs():
	'''
	==============================================================================================================
	When called, function collects input arguments from the command line. Outputs variables to be used by main()
	==============================================================================================================
	'''

	# Defining arguments
=======
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
from matplotlib.colors import LogNorm
from matplotlib import rc
matplotlib.rcParams['text.usetex'] = True		# Enable LaTeX in annotations. 
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{siunitx} \sisetup{detect-all} \usepackage{helvet} \usepackage{sansmath} \sansmath' # Sans-serif math
from scipy import ndimage, misc
import cv2
import sys
import argparse


def CollctUserArgs():
	# Arguments
>>>>>>> ad95adbefe0997511532561a48002d6b2517df26
	parser = argparse.ArgumentParser()

	parser.add_argument('-v', "--verbose", help="Increase output verbosity",
	                    action="store_true")

<<<<<<< HEAD
	parser.add_argument('-p', "--pdb", help="PDB file name. Include path if necessary.",
						type=str)

	parser.add_argument('-m', "--models", help="Desired model number (default = use all models). Model number corresponds to order in PDB file.",
						type=int)

	parser.add_argument('-c', "--chains", help="Desired chain number (default = use all chains). Chain number corresponds to order in PDB file.",
						type=int)

	parser.add_argument('-d', "--out_dir", help="Out directory. Must be available before-hand.",
                    	type=str)

	parser.add_argument('-t', "--plot_type", help="Type of angles plotted Ramachandran diagram. Refer to README.md for options and details.",
	                    type=int)

	parser.add_argument('-s', "--save_csv", help="Saves calculated dihedral angles in a separate CSV file.",
	                    action="store_true")

	args = parser.parse_args()

	# Analysing arguments 
	if not args.models:			# Iterate models?
=======
	parser.add_argument('-p', "--pdb", help="PDB file name. Include path if necessary.")

	parser.add_argument('-m', "--models", help="Desired model number. Model number corresponds to order in PDB file. Must set --iter_models to False first.",
						type=int)

	parser.add_argument('-c', "--chains", help="Desired chain number. Chain number corresponds to order in PDB file. Must set --iter_chains to False first.")

	parser.add_argument('-d', "--out_dir", help="Out directory. Must be available beforehand.",
                    	type=str)

	parser.add_argument('-t', "--plot_type", help="Type of reference angles used in Ramachandran plot(s). Refer to documentation for options and details.",
	                    type=int)

	args = parser.parse_args()

	# Analysing arguments 

	# Iterate models?
	if not args.models:
>>>>>>> ad95adbefe0997511532561a48002d6b2517df26
		model_num = 0
		itmod = True
	elif isinstance(args.models, int):
		model_num = args.models
		itmod = False
	else:
		print('Invalid model number. ')

<<<<<<< HEAD
	if not args.chains:			# Iterate chains?
=======
	# Iterate chains?
	if not args.chains:
>>>>>>> ad95adbefe0997511532561a48002d6b2517df26
		chain_num = 0
		itchain = True
	elif isinstance(args.chains, int):
		chain_num = args.chains
		itchain = False
	else:
		print('Invalid chain number.')

<<<<<<< HEAD
	if not args.out_dir:		# Handling the out directory
=======
	# Handling the out directory
	if not args.out_dir:
>>>>>>> ad95adbefe0997511532561a48002d6b2517df26
		out_dir = './'
	elif isinstance(args.out_dir, str):
		out_dir = args.out_dir
	else:
		print('Enter valid directory to write Ramachandran plot.')

<<<<<<< HEAD
	if args.plot_type is None:
		plot_type = 0
	else:
		plot_type = args.plot_type
		print('No plot type given. Defaulting to plot all.')

	return args.pdb, itmod, model_num, itchain, chain_num, plot_type, out_dir, args.verbose, args.save_csv



def VerboseStatement(verb_boolean, statement):
	'''
	=============================================================================
	If first argument is true, function prints a given statement to command line.
	=============================================================================
	'''

	if verb_boolean: 		# Parsed from argparser
		print(statement)
		# print('\n')
=======
	return args.pdb, itmod, model_num, itchain, chain_num, args.plot_type, out_dir, args.verbose


def VerboseStatement(verb_boolean, statement):
	if verb_boolean:
		print(statement)
>>>>>>> ad95adbefe0997511532561a48002d6b2517df26
	else:
		pass