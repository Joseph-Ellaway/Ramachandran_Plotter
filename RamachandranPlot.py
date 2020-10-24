
''' ================================================================================================================

	Plots the frequency density of the dihedral (torsion) angles from the Top8000 peptide database as a background 
	histogram, along with two contour lines. 

	The user's input peptide's dihedral angles are then added ontop of this background as a scatter plot. 

	Colours and recommended parameters can be easily adjusted towards the end of the script. 

	================================================================================================================ '''


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

background_colour = 'Blues'



########################################################
#						FUNCTIONS 					   #
########################################################

def PhiPsiDataFramer(dataframe_dict):						# dtypes : dict of Pandas DataFrames
	'''
	=======================================================
	Takes a dictionary of Pandas DataFrames and extracrs phi and psi angles corresponding to the bond type in each DataFrame.
	=======================================================
	'''
	phi_psi_dict = {}
	for dataframe_key, dataframe_value in dataframe_dict.items():
		exec('phi_psi_dict.update({"%s_angles" : [dataframe_value["phi"], dataframe_value["psi"]]})' % (dataframe_key,))
	return phi_psi_dict



def AxesRemover(mpl_axis):
	'''
	=======================================================
	Removes axes and tick marks from a given Matplotlib axis. 
	=======================================================
	'''
	mpl_axis.set_frame_on(False)
	mpl_axis.get_xmpl_axis().tick_bottom()
	mpl_axis.get_ympl_axis().set_visible(False)
	mpl_axis.get_xmpl_axis().set_visible(False)



def SaveAndCloseFigure(out_file_name, resolution):			# dtypes : string, int 
	'''
	=======================================================
	Saves the plot area (defined before function called) as a PNG image to given file name and resolution. 
	=======================================================
	'''
	plt.savefig(out_file_name, dpi=resolution, bbox_inches=0, pad_inches=None)
	plt.close()



def PhiPsiPlotter(phi_angles, psi_angles, name, figsize):	# dtypes : array-like, 
	'''
	=======================================================
	Creates a 2D histogram (Ramachandran plot) of out of phi and psi angles from the Top8000 peptide DB.
	=======================================================
	'''
	fig, ax = plt.subplots(1,1, figsize=figsize, tight_layout=True)
	fig.patch.set_visible(False)
	kwargs = {
		'bins' : 140, 
		'cmap' : background_colour, 
		'norm' : colors.PowerNorm(0.1),
		'alpha' : 1}
	# Reference https://matplotlib.org/3.2.1/gallery/statistics/hist.html for help with making 2D histogram.
	ax.set_axis_off()
	ax.hist2d(phi_angles, psi_angles, **kwargs)
	out_file_name = str(name[:-7] + 'RamachangranPlot.png')
	SaveAndCloseFigure(out_file_name, 80)
	return out_file_name



def Smoother(out_file_name, figsize):
	'''
	=======================================================
	Reads a given pixelated image (passed into function as the file name as a string) and 'smoothes' it out to create an unpixelated version.
	The unpixelated ('smoothed') version is then saved to ~/your_current_dir, overwriting the pixelated input image. 
	=======================================================
	'''
	fig, ax = plt.subplots(1,1, figsize=figsize, tight_layout=True)
	fig.patch.set_visible(False)
	kwargs = {
		'cmap' : str(background_colour + '_r'), 
		'alpha' : 1}
	ax.set_axis_off()
	rama_plot = cv2.imread(out_file_name, 0)
	blurred_rama_plot = ndimage.gaussian_filter(rama_plot, sigma=0.3)
	smoothed_rama_plot = ndimage.percentile_filter(blurred_rama_plot, percentile=90, size=20)
	ax.imshow(smoothed_rama_plot, **kwargs)
	bkgd_resolution = 96									# Adjust integer to change background resolution (higher = better quality but slower run time)
	SaveAndCloseFigure(out_file_name, bkgd_resolution)
	plt.close()
	





########################################################
#					IMPORTING DATA 					   #
########################################################

# CONDITIONAL CHECK TO ENSURE USER SELECTS THE DESIRED RAMACHANDRAN PLOT
if len(sys.argv) < 3:
	print('\n \n Enter valid PDB file name and Ramachandran plot code. \n Available options are: \n \
#################################### \n \
	0 = All angles \n \
	1 = General angles \n \
	2 = Glycine angles \n \
	3 = Trans-proline angles only \n \
	4 = Cis-proline angles only \n \
	5 = Proline (cis/trans) angles only \n \
	6 = Pre-proline angles only \n \
	7 = Ile-Val angles only \n \
#################################### \n')
	print('\n \n python3 RamachangranPlot.py <pdb-file-name> <plot-code> \n \n e.g. \n python3 RamachangranPlot.py pep1.pdb 2 \n \n')
	sys.exit()


# TOP8000 GOLD STANDARD DATA
top8000_df = pd.read_csv('Top8000_DihedralAngles.csv')		# Processed dihedral angle data from Top8000 peptide data set. 

dataframe_dict = {			# Top8000 peptide data classified into keys based on bond type. Only remove entries if you are experiencing memory issues. 
	'All' : top8000_df,
	'General' : top8000_df.loc[top8000_df['type']=='General'], 
	'Glycine' : top8000_df.loc[top8000_df['type']=='Glycine'], 
	'Trans-proline' : top8000_df.loc[top8000_df['type']=='Trans-proline'],
	'Cis-proline' : top8000_df.loc[top8000_df['type']=='Cis-proline'],
	'Proline' : top8000_df.loc[top8000_df['type']=='Cis-proline'].append(top8000_df.loc[top8000_df['type']=='Trans-proline']), 
	'Pre-proline' : top8000_df.loc[top8000_df['type']=='Pre-proline'], 
	'Ile-Val' : top8000_df.loc[top8000_df['type']=='Isoleucine']
	}

top8000_phipsi_dict = PhiPsiDataFramer(dataframe_dict)				# Extracting phi-psi angle data. 


# USER'S INPUT DIHEDRAL ANGLE DATA
dihedral_angle_data = pd.read_csv(str(sys.argv[1][:-4] + '_DihedralAngles.csv'))		# Name of .csv with user's dihedral angle information. 

favoured_dihedral_angles = dihedral_angle_data[(dihedral_angle_data['evaluation'] == 'Favored') | (dihedral_angle_data['evaluation'] == 'Allowed')]
outlier_dihedral_angles = dihedral_angle_data[(dihedral_angle_data['evaluation'] == 'OUTLIER')]

user_favoured_dataframe_dict = {			# Top8000 peptide data classified into keys based on bond type. Only remove entries if you are experiencing memory issues. 
	'All' : favoured_dihedral_angles,
	'General' : favoured_dihedral_angles.loc[favoured_dihedral_angles['type']=='General'], 
	'Glycine' : favoured_dihedral_angles.loc[favoured_dihedral_angles['type']=='Glycine'], 
	'Trans-proline' : favoured_dihedral_angles.loc[favoured_dihedral_angles['type']=='Trans-proline'],
	'Cis-proline' : favoured_dihedral_angles.loc[favoured_dihedral_angles['type']=='Cis-proline'],
	'Proline' : favoured_dihedral_angles.loc[favoured_dihedral_angles['type']=='Cis-proline'].append(favoured_dihedral_angles.loc[favoured_dihedral_angles['type']=='Trans-proline']), 
	'Pre-proline' : favoured_dihedral_angles.loc[favoured_dihedral_angles['type']=='Pre-proline'], 
	'Ile-Val' : favoured_dihedral_angles.loc[favoured_dihedral_angles['type']=='Isoleucine']
	}

user_outlier_dataframe_dict = {			# Top8000 peptide data classified into keys based on bond type. Only remove entries if you are experiencing memory issues. 
	'All' : outlier_dihedral_angles,
	'General' : outlier_dihedral_angles.loc[outlier_dihedral_angles['type']=='General'], 
	'Glycine' : outlier_dihedral_angles.loc[outlier_dihedral_angles['type']=='Glycine'], 
	'Trans-proline' : outlier_dihedral_angles.loc[outlier_dihedral_angles['type']=='Trans-proline'],
	'Cis-proline' : outlier_dihedral_angles.loc[outlier_dihedral_angles['type']=='Cis-proline'],
	'Proline' : outlier_dihedral_angles.loc[outlier_dihedral_angles['type']=='Cis-proline'].append(outlier_dihedral_angles.loc[outlier_dihedral_angles['type']=='Trans-proline']), 
	'Pre-proline' : outlier_dihedral_angles.loc[outlier_dihedral_angles['type']=='Pre-proline'], 
	'Ile-Val' : outlier_dihedral_angles.loc[outlier_dihedral_angles['type']=='Isoleucine']
	}







########################################################
#			GENERATING FAVOURED REGION(S)			   #
########################################################

figure_size_background=(10,10)								# Does not change output image size. 
options = ['All', 'General', 'Glycine', 'Trans-proline', 'Cis-proline', 'Proline', 'Pre-proline', 'Ile-Val']

chosen_plot = int(sys.argv[2])								# User input determines what Ramachandran background is plotted. 
angle_type = str(options[chosen_plot] + '_angles')
rama_plot_name = PhiPsiPlotter(top8000_phipsi_dict.get(angle_type)[0], top8000_phipsi_dict.get(angle_type)[1], angle_type, figure_size_background)
Smoother(rama_plot_name, figure_size_background)


# EXTRACTING DIHEDRAL ANGLES FROM USER'S DATA
favoured_phipsi_dict = PhiPsiDataFramer(user_favoured_dataframe_dict)
outlier_phipsi_dict = PhiPsiDataFramer(user_outlier_dataframe_dict)

user_chosen_favoured_dihedral_angles = favoured_phipsi_dict.get(angle_type)
user_chosen_outlier_dihedral_angles = outlier_phipsi_dict.get(angle_type)






########################################################
#					PLOTTING DATA					   #
########################################################

# AESTHETIC CHOICES
plt.style.use('seaborn-poster')
zero_lines_kwargs = {
	'colors' : ['Grey'], 
	'alpha' : 0.4, 
	'zorder' : 3, 
	'linewidths' : [1]
	}
# Recommended adjustable parameters. 
figure_size = (8,8)						# Output figure size in inches.
contour_line_color_inner = '#DFF8FB'	# Inner contour line colour.
contour_level_inner = 96				# Percentile of dihedral angles for inner contour lines (e.g. contour_level=96 means the area bounded by the contour line represents the range of angles in which 96% of all dihedral from the Top800 peptide DB fall within).
contour_level_outer = 15				# Percentile of dihedral angles for outer contour lines
contour_line_color_outer = '#045E93'	# Colour of outer contour lines
contour_line_alpha_outer = 0.2			# Opacity of outer contour lines
out_resolution = 200					# Output figure resolution.
data_point_colour = '#D4AB2D'			# Colour of data points for each Phi-Psi dihedral angle pair. 
data_point_edge_colour = '#3c3c3c'		# Colour of data point's border.
outlier_colour = 'Red'					# Colour of data point outliers. 

fig, ax = plt.subplots(1,1, figsize=figure_size, tight_layout=True)		# Defining plot area. 


# ADDING COUNTOURS - Comment this section out to remove contour lines from plot area.
if chosen_plot == 4:																	# Reduced number of bins for Cis-Pro Ramachandran plot.
	counts, discard1, discard2, discard3 = plt.hist2d(top8000_phipsi_dict.get(angle_type)[0],top8000_phipsi_dict.get(angle_type)[1], bins=45, norm=LogNorm(), alpha=0)
	ax.contour(counts.transpose(), extent=[-180, 180, -180, 180], levels=[contour_level_inner], linewidths=1, colors=[contour_line_color_inner], zorder=2)

	counts_v2, discard1, discard2, discard3 = plt.hist2d(top8000_phipsi_dict.get(angle_type)[0],top8000_phipsi_dict.get(angle_type)[1], bins=35, norm=LogNorm(), alpha=0)
	ax.contour(counts_v2.transpose(), extent=[-180, 180, -180, 180], levels=[contour_level_outer], linewidths=1, colors=[contour_line_color_outer], zorder=2, alpha=contour_line_alpha_outer)

else:
	counts, discard1, discard2, discard3 = plt.hist2d(top8000_phipsi_dict.get(angle_type)[0],top8000_phipsi_dict.get(angle_type)[1], bins=90, norm=LogNorm(), alpha=0)
	ax.contour(counts.transpose(), extent=[-180, 180, -180, 180], levels=[contour_level_inner], linewidths=1, colors=[contour_line_color_inner], zorder=2)

	counts_v2, discard1, discard2, discard3 = plt.hist2d(top8000_phipsi_dict.get(angle_type)[0],top8000_phipsi_dict.get(angle_type)[1], bins=35, norm=LogNorm(), alpha=0)
	ax.contour(counts_v2.transpose(), extent=[-180, 180, -180, 180], levels=[contour_level_outer], linewidths=1, colors=[contour_line_color_outer], zorder=2, alpha=contour_line_alpha_outer)

# ADDING FAVOURED RAMACHANDRAN REGION IMAGE TO BACKGROUND 
ax.imshow(plt.imread(rama_plot_name), extent=[-195, 195, -195, 195], zorder=1)


# ADDING GRIDLINES
ax.grid(alpha=0.4, linestyle='--', linewidth=1, color='Grey', zorder=3)
ax.hlines(0, -180, 180, **zero_lines_kwargs)
ax.vlines(0, -180, 180, **zero_lines_kwargs)

# PLOTTING USER'S DIHEDRAL ANGLE DATA
ax.scatter(user_chosen_outlier_dihedral_angles[0], user_chosen_outlier_dihedral_angles[1], s=15, color=outlier_colour, zorder=4, linewidths=0.5, edgecolor=data_point_edge_colour)
ax.scatter(user_chosen_favoured_dihedral_angles[0], user_chosen_favoured_dihedral_angles[1], s=15, color=data_point_colour, zorder=4, linewidths=0.5, edgecolor=data_point_edge_colour)


# AXES AESTHETICS/FEATURES
ax.set_xlim((-180, 180))
ax.set_ylim((-180, 180))
ax.set_xlabel(r'$\phi \ (\si{\degree})$')
ax.set_ylabel(r'$\psi \ (\si{\degree})$')
ax_linewidth = 2
ax.spines['left'].set_linewidth(ax_linewidth)
ax.spines['bottom'].set_linewidth(ax_linewidth)


# SAVING RAMACHANDRAN PLOT AS PNG IMAGE
plt.savefig(rama_plot_name, dpi=out_resolution)	# dpi = resolution. Adjust to suite your prefernces
print('\n Done. \n Ramachangran plot saved to', rama_plot_name, '\n \n')
