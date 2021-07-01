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
	# parser.add_argument('-v', "--verbose", help="Increase output verbosity",
	#                     action="store_true")
	parser.add_argument("--pdb", help="PDB file name. Include path if necessary.")
	# parser.add_argument("--iter_models", help="If PDB file contains multiple models, calculate Ramachandran angles of all models if True (default). Set to False if a single/selection of models is required.",
	#                     action="store_true")
	parser.add_argument("--models", help="Desired model number. Model number corresponds to order in PDB file. Must set --iter_models to False first.",
	                    action="store_true")
	# parser.add_argument("--iter_chains", help="If model in PDB file contains multiple chains, calculate Ramachandran angles for all chains if True (default). Set to False if a single/selection of chains is required.",
	#                     action="store_true")
	parser.add_argument("--chains", help="Desired chain number. Chain number corresponds to order in PDB file. Must set --iter_chains to False first.",
	                    action="store_true")
	parser.add_argument("--out_dir", help="Out directory. Must be available beforehand.",
	                    action="store_true")
	parser.add_argument("--plot_type", help="Type of reference angles used in Ramachandran plot(s). Refer to documentation for options and details.",
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
		itchn = True
	elif isinstance(args.chains, int):
		chain_num = args.chains
		itchn = False
	else:
		print('Invalid chain number.')


	return args.pdb, itmod, model_num, itchain, chain_num, plot_type