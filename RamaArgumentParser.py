""" 
	====================================================================================
	If running script from the command line, functions here are called to parsed user's 
	arguments into the main() function in RamachandranPlotter.py.
	
	Version 2.0.1:
	 - Relies on the easily accessible Biopython package, rather than Phenix as in 
	   versions <2.0
	 - User arguments can be now easily parsed in from the command line (as facilitated 
	   by functions here)
	 - If required, the script could be implemented into existing protein analysis 
	   pipelines by importing this function ( main() ).

	Author information:
	 - Joseph I. J. Ellaway
	 - josephellaway@gmail.com
	 - https://github.com/Joseph-Ellaway
	====================================================================================
"""

import argparse



def CollctUserArgs():
	"""
	====================================================================================
	When called, function collects input arguments from the command line. Outputs 
	variables to be used by main()
	====================================================================================
	"""

	# Defining arguments
	parser = argparse.ArgumentParser()

	parser.add_argument("-v", "--verbose", help="Increase output verbosity",
	                    action="store_true")

	parser.add_argument("-p", "--pdb", 
						help="PDB file name: <filename.pdb>. Include path if necessary.",
						type=str)

	parser.add_argument("-m", "--models", 
						help="Desired model number (default: all models). Model number corresponds to order in PDB file.",
						type=int)

	parser.add_argument("-c", "--chains", 
						help="Desired chain number (default: all chains). Chain number corresponds to order in PDB file.",
						type=int)

	parser.add_argument("-d", "--out_dir", 
						help="Out directory. Must be available before-hand.",
                    	type=str)

	parser.add_argument("-t", "--plot_type", 
						help="Type of angles plotted on Ramachandran diagram. Refer to README.md for options and details.",
	                    type=int)

	parser.add_argument("-f", "--file_type", 
						help="File type for output plot. Options: PNG (default, 96 dpi), PDF, SVG, EPS and PS.",
						type=str)

	parser.add_argument("-s", "--save_csv", 
						help="Save calculated dihedral angles in separate CSV.",
	                    action="store_true")

	args = parser.parse_args()

	# Analysing arguments 
	if not args.models:			# Iterate models?
		model_num = 0
		itmod = True
	elif isinstance(args.models, int):
		model_num = args.models
		itmod = False
	else:
		print("Invalid model number. ")

	if not args.chains:			# Iterate chains?
		chain_num = 0
		itchain = True
	elif isinstance(args.chains, int):
		chain_num = args.chains
		itchain = False
	else:
		print("Invalid chain number.")

	if not args.out_dir:		# Handling the out directory
		out_dir = "./"
	elif isinstance(args.out_dir, str):
		out_dir = args.out_dir
	else:
		print("Enter valid directory to write Ramachandran plot.")

	if args.plot_type is None:
		plot_type = 0
	elif args.plot_type in [0, 1, 2, 3, 4, 5]:
		plot_type = args.plot_type
	else:
		print("Invalid plot type given. Give integer value between 0-5 to compute dihedral angles.")
		print("	E.g. 	--plot_type <int>")
		print("Options:")
		print("	0 : All \n \
	1 : General (All residues bar Gly, Pro, Ile, Val and pre-Pro) \n \
	2 : Glycine \n \
	3 : Proline (cis and trans) \n \
	4 : Pre-proline (residues preceeding a proline) \n \
	5 : Ile or Val")
		exit()

	if not args.file_type:
		file_type = "png"
	else:
		# convert file type to lower case
		file_type = args.file_type.lower()

	return args.pdb, itmod, model_num, itchain, chain_num, plot_type, out_dir, args.verbose, args.save_csv, file_type



def VerboseStatement(verb_boolean, statement):
	"""
	=============================================================================
	If first argument is true, function prints a given statement to command line.
	=============================================================================
	"""

	if verb_boolean: 		# Parsed from argparser
		print(statement)
		# print("\n")
	else:
		pass