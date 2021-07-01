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
	parser = argparse.ArgumentParser()

	parser.add_argument('-v', "--verbose", help="Increase output verbosity",
	                    action="store_true")

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
		model_num = 0
		itmod = True
	elif isinstance(args.models, int):
		model_num = args.models
		itmod = False
	else:
		print('Invalid model number. ')

	# Iterate chains?
	if not args.chains:
		chain_num = 0
		itchain = True
	elif isinstance(args.chains, int):
		chain_num = args.chains
		itchain = False
	else:
		print('Invalid chain number.')

	# Handling the out directory
	if not args.out_dir:
		out_dir = './'
	elif isinstance(args.out_dir, str):
		out_dir = args.out_dir
	else:
		print('Enter valid directory to write Ramachandran plot.')

	return args.pdb, itmod, model_num, itchain, chain_num, args.plot_type, out_dir, args.verbose


def VerboseStatement(verb_boolean, statement):
	if verb_boolean:
		print(statement)
	else:
		pass