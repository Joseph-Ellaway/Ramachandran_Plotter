
''' ================================================================================================================

	Plots the frequency density of the dihedral (torsion) angles from the Top8000 peptide database as a background 
	histogram, along with two contour lines. 

	The user's input peptide's dihedral angles are then added ontop of this background as a scatter plot. 

	Colours and recommended parameters can be easily adjusted towards the end of the script. 

	================================================================================================================ '''

import pandas as pd
import matplotlib.pyplot as plt

from DihedralCalculator import *
from PlotterFunctions import *
from RamaArgumentParser import *



# Arguments
pdb, itmod, model_num, itchain, chain_num, plot_type, out_dir, verb = CollctUserArgs()



########################################################
#				IMPORTING USER DATA					   #
########################################################

VerboseStatement(verb, str('Importing ' + pdb + '\n') )

userpdb_df = ExtractDihedrals(pdb_file_name=pdb, iter_models=itmod, model_number=model_num, iter_chains=itchain, chain_id=chain_num)
userpdb_df = userpdb_df.dropna()

VerboseStatement(verb, 'Dihedral angles calculated')




# ########################################################
# #				IMPORTING REFERENCE DATA			   #
# ########################################################

VerboseStatement(verb, 'Importing Top8000 library \n')

# TOP8000 GOLD STANDARD DATA
top8000_df = pd.read_csv('Top8000_DihedralAngles.csv')		# Processed dihedral angle data from Top8000 peptide data set. 

# Uses Top8000 DB -- generating background
options = ['All', 'General', 'Glycine', 'Proline', 'Pre-proline']
plot_type = options[int(plot_type)]								# User input determines what Ramachandran background is plotted. 
plot_name = str(plot_type + 'RamachangranPlot.png')

VerboseStatement(verb, 'Generating background of favoured regions \n')

MakeBackground(top8000_df, plot_type, plot_name)



########################################################
#					PLOTTING DATA					   #
########################################################
VerboseStatement(verb, 'Plotting Ramachandran diagram \n')


plt.style.use('seaborn-poster')


# Recommended adjustable parameters. 
figure_size = (8,8)						# Output figure size in inches.
contour_level_inner = 96				# Percentile of dihedral angles for inner contour lines (e.g. contour_level=96 means the area bounded by the contour line represents the range of angles in which 96% of all dihedral from the Top800 peptide DB fall within).
contour_level_outer = 15				# Percentile of dihedral angles for outer contour lines
contour_line_color_inner = '#DFF8FB'	# Inner contour line colour.
contour_line_color_outer = '#045E93'	# Colour of outer contour lines
out_resolution = 200					# Output figure resolution.
data_point_colour = '#D4AB2D'			# Colour of data points for each Phi-Psi dihedral angle pair. 
data_point_edge_colour = '#3c3c3c'		# Colour of data point's border.

fig, ax = plt.subplots(1,1, figsize=figure_size, tight_layout=True)		# Defining plot area. 


# # ADDING COUNTOURS - Comment this section out to remove contour lines from plot area.
AddContour(ax, top8000_df, contour_level=contour_level_inner, line_colour=contour_line_color_inner)
AddContour(ax, top8000_df, contour_level=contour_level_outer, line_colour=contour_line_color_outer, contour_alpha=0.3)


# ADDING FAVOURED RAMACHANDRAN REGION IMAGE TO BACKGROUND 
ax.imshow(plt.imread(plot_name), extent=[-195, 195, -195, 195], zorder=1)


# ADDING GRIDLINES
AddGridLines(ax)

# PLOTTING USER'S DIHEDRAL ANGLE DATA
ax.scatter(userpdb_df['phi'], userpdb_df['psi'], s=15, color=data_point_colour, zorder=4, linewidths=0.5, edgecolor=data_point_edge_colour)


# AXES AESTHETICS/FEATURES
FormatAxis(ax)


# SAVING RAMACHANDRAN PLOT AS PNG IMAGE
SaveAndCloseFigure(plot_name, out_resolution)
plt.show()
print('\n Done. \n Ramachangran plot saved to', plot_name, '\n \n')
