import argparse
import logging

logger = logging.getLogger(__name__)

PLOT_TYPES_ERROR_MSG = """Invalid plot type given.

Give integer value between 0-5 to compute dihedral angles.

Options:
	"all" 			: All
	"general" 		: General (All residues bar Gly, Pro, Ile, Val and pre-Pro)
	"glycine" 		: Glycine
	"proline" 		: Pro
	"pre-proline" 	: Pre-proline (residues preceeding a proline)
	"ile-val" 		: Ile or Val
"""

def parse_user_args():
	"""
	When called, function collects input arguments from the command line. Outputs
	variables to be used by main()
	"""

	# Defining arguments
	parser = argparse.ArgumentParser()

	parser.add_argument(
		"-v",
		"--verbose",
		help="Increase output verbosity",
	    action="store_true"
	)

	parser.add_argument(
		"-i",
		"--input-file",
		help="""
			Path(s) to mmCIF file(s). Separate multiple files with a comma.
			E.g. file1.cif,file2.cif.

			To specify chains, use multiple -i flags for each file, followed by the
			label_asym_id chain ID. E.g. -i file1.cif:A,file2.cif:B,C
			""",
		type=str,
        nargs="+",
        action="append",
        required=True,
	)

	parser.add_argument(
		"-o",
		"--out-plot-file",
		help="Directory and filename of output plot.",
        type=str,
		required=True
	)

	parser.add_argument(
		"-p",
		"--plot-type",
		default="all",
		help="Type of angles plotted on Ramachandran diagram. Refer to README.md for options and details.",
	    type=str
	)

	parser.add_argument(
		"-f",
		"--file-name",
		help="File name for output plot. Defaults to 'ramachandran_plot'.",
		type=str,
		default="ramachandran_plot"
	)

	parser.add_argument(
		"-t",
		"--file-type",
		help="File type for output plot. Options: PNG (default, 96 dpi), PDF, SVG, EPS and PS.",
		default="png",
		type=str
	)

	parser.add_argument(
		"-s",
		"--save-csv-file",
		help="Path to save calculated dihedral angles in separate CSV.",
		type=str,
	)

	parser.add_argument(
		"-r",
		"--remove-outliers",
		help="Remove outliers from the plot. Determined by representative angles.",
		action="store_true",
		default=False
	)

	parser.add_argument(
		"-u",
		"--use-struct-asym-ids",
		help="Use structure's asym IDs as chain IDs, rather than author-defined chain IDs.",
		action="store_true",
		default=False
	)

	parser.add_argument(
		"-e",
		"--existing-rama-background",
		help="Path to existing Ramachandran distribution background file.",
		type=str,
		default=None
	)

	parser.add_argument(
		"-x",
		"--exclude-residues",
		help="Exclude additional residues from the plot.",
		type=list,
		default=[]
	)

	parser.add_argument(
		"-c",
		"--canonical-only",
		default=False,
		help="Only plot canonical residues.",
		action="store_true"
	)

	parser.add_argument(
		"-n",
		"--exclude-pre-proline",
		help="Exclude residues processing pre-proline residues.",
		action="store_true",
		default=False
	)

	args = parser.parse_args()

	# # Analysing arguments
	# if not args.models:			# Iterate models?
	# 	model_num = 0
	# 	itmod = True
	# elif isinstance(args.models, int):
	# 	model_num = args.models
	# 	itmod = False
	# else:
	# 	print("Invalid model number. ")

	# if not args.chains:			# Iterate chains?
	# 	chain_num = 0
	# 	itchain = True

	# elif isinstance(args.chains, int):
	# 	chain_num = args.chains
	# 	itchain = False
	# else:
	# 	print("Invalid chain number.")

	if args.plot_type not in (
		"all","general", "glycine","proline","pre-proline","ile-val"
	):

		logger.error(PLOT_TYPES_ERROR_MSG)

		raise ValueError("Invalid plot type given. ")

	return args
