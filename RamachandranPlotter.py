
''' ================================================================================================================

	Plots the frequency density of the dihedral (torsion/Ramachandran) angles from the Top8000 peptide database as a 
	smoothed background histogram, along with two contour lines. 

	The user's input peptide's dihedral angles are then added ontop of this background as a scatter plot. 

	Colours and recommended parameters can be easily adjusted towards the end of the script. 
	
	Version 2.0.1:
	 - Relies on the easily accessible Biopython package, rather than Phenix as in versions <2.0
	 - User arguments can be now easily parsed in from the command line
	 - If required, the script could be implemented into existing protein analysis pipelines by importing this 
	 function ( main() )

	Author information:
	 - Joseph I. J. Ellaway
	 - josephellaway@gmail.com
	 - https://github.com/Joseph-Ellaway

	================================================================================================================ '''

# Base functions
import pandas as pd
import matplotlib.pyplot as plt

# Package functions
from DihedralCalculator import *
from PlotterFunctions import *
from RamaArgumentParser import *

# Main function
def main(pdb, itmod, model_num, itchain, chain_num, plot_type, out_dir, verb, save):

	########################################################
	#				IMPORTING USER DATA					   #

	VerboseStatement(verb, str('Importing ' + str(pdb)) )

	userpdb_df = ExtractDihedrals(pdb_file_name=pdb, iter_models=itmod, model_number=model_num, iter_chains=itchain, chain_id=chain_num)
	userpdb_df = userpdb_df.dropna()   # Removes invalid dihedral angles and angles from ligands/non-canonical amino acid residues

	# Selecting user's desired Ramachandran plot
	options = ['All', 'General', 'Glycine', 'Proline', 'Pre-proline', 'Ile-Val'] # Available Ramachandran plots
	plot_type = options[int(plot_type)]								  # User input determines what Ramachandran background is plotted. 
	plot_name = str(out_dir + plot_type + 'RamachandranPlot')	  # Out file name



	########################################################
	#				IMPORTING REFERENCE DATA			   #

	VerboseStatement(verb, 'Importing Top8000 library')

	# Top8000 peptide dataset. Pre-analysed
	top8000_df = pd.read_csv('Top8000_DihedralAngles.csv.gz', compression='gzip') 



	########################################################
	#				SELECTING RESIDUE TYPE DATA			   #

	if plot_type == 'All':
		pass

	else:
		userpdb_df = userpdb_df.loc[userpdb_df['type'] == plot_type]



	########################################################
	#				SAVING USER DATA (optional)			   #

	if save:
		csv_file_name = str(plot_name + '.csv')
		VerboseStatement(verb, str('Saving CSV as: ' + csv_file_name))
		userpdb_df.to_csv(csv_file_name, index=False)

	else:
		pass

	VerboseStatement(verb, 'Dihedral angles calculated')



	########################################################
	#					PLOTTING DATA					   #

	# Recommended adjustable parameters. 
	figure_size = (5,5)						# Output figure size in inches.
	contour_level_inner = 96				# Percentile of dihedral angles for inner contour lines (e.g. contour_level=96 means the area bounded by the contour line represents the range of angles in which 96% of all dihedral from the Top800 peptide DB fall within).
	contour_level_outer = 15				# Percentile of dihedral angles for outer contour lines
	contour_line_color_inner = '#DFF8FB'	# Inner contour line colour.
	contour_line_color_outer = '#045E93'	# Colour of outer contour lines
	out_resolution = 200					# Output figure resolution. Not required if saving file as PDF
	data_point_colour = '#D4AB2D'			# Colour of data points for each Phi-Psi dihedral angle pair. 
	data_point_edge_colour = '#3c3c3c'		# Colour of data point's border.
	background_colour = 'Blues'				# Colour map of background plot. Refer to https://matplotlib.org/stable/tutorials/colors/colormaps.html for colormap options

	# Plotting background
	VerboseStatement(verb, 'Generating background of favoured regions')

	MakeBackground(top8000_df, plot_type, plot_name, background_colour) 	# Genertating background: region of favoured dihedral angles

	# Plotting user's PDB dihedral angles
	VerboseStatement(verb, 'Plotting Ramachandran diagram')

	plt.style.use('seaborn-poster')

	fig, ax = plt.subplots(1,1, figsize=figure_size, tight_layout=True)		# Defining plot area. 

	# ADDING COUNTOURS - Comment this section out to remove contour lines from plot area.
	AddContour(ax, top8000_df, contour_level=contour_level_inner, line_colour=contour_line_color_inner)
	AddContour(ax, top8000_df, contour_level=contour_level_outer, line_colour=contour_line_color_outer, contour_alpha=0.3)

	# ADDING FAVOURED RAMACHANDRAN REGION IMAGE TO BACKGROUND 
	ax.imshow(plt.imread(str(plot_name + '.png')), extent=[-195, 195, -195, 195], zorder=1)

	# ADDING GRIDLINES
	AddGridLines(ax)

	# PLOTTING USER'S DIHEDRAL ANGLE DATA
	ax.scatter(userpdb_df['phi'], userpdb_df['psi'], s=15, color=data_point_colour, zorder=4, linewidths=0.5, edgecolor=data_point_edge_colour)

	# AXES AESTHETICS/FEATURES
	FormatAxis(ax)

	# SAVING RAMACHANDRAN PLOT AS PNG IMAGE
	VerboseStatement(verb, 'Saving plot')

	# ... as PNG
	plt.savefig(str(plot_name + '.png'), dpi=out_resolution, bbox_inches=0, pad_inches=None)

	# ... as PDF
	plt.savefig(str(plot_name + '.pdf'), bbox_inches=0, pad_inches=None)
	plt.close()

	# plt.show()

	print('Done. \n Ramachandran plot saved to', str(plot_name + '.png'))





########################################################
#					RUNNING SCRIPT					   #

if __name__ == "__main__":

	# Loading user's input arguments
	pdb, itmod, model_num, itchain, chain_num, plot_type, out_dir, verb, save = CollctUserArgs()

	main(pdb, itmod, model_num, itchain, chain_num, plot_type, out_dir, verb, save)

else:
	pass