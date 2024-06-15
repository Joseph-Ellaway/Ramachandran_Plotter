
"""
    ====================================================================================
    Plots the frequency density of the dihedral (torsion/Ramachandran) angles from the
    Top8000 peptide database as a smoothed background histogram, along with two contour
    lines.

    The user's input peptide's dihedral angles are then added ontop of this background
    as a scatter plot.

    Colours and recommended parameters can be easily adjusted towards the end of the
    script.

    Version 2.0.1:
     - Relies on the easily accessible Biopython package, rather than Phenix as in
       versions <2.0
     - User arguments can be now easily parsed in from the command line
     - If required, the script could be implemented into existing protein analysis
       pipelines by importing this function ( main() )

    Author information:
     - Joseph I. J. Ellaway
     - josephellaway@gmail.com
     - https://github.com/Joseph-Ellaway
    ====================================================================================
"""

import os
import matplotlib.pyplot as plt
import pandas as pd
import logging
import json

# Package functions
from .rama_calc_utils import *
from .plot_utils import *
from .argument_parser import *

# Logging
logger = logging.getLogger(__name__)



def calculate_rama_angles(pdb_file, chain_id=None):
    """
    Calculate the phi and psi dihedral angles of a PDB file
    """
    # Extract dihedral angles
    dihedral_angles = extract_dihedrals(
        pdb_file_name=pdb_file,
        chain_id=chain_id
    )

    # Remove invalid dihedral angles/angles from ligands or non-canonical residues
    dihedral_angles = dihedral_angles.dropna()

    return dihedral_angles




class RamachandranPlotter:
    """
    Object to parse, calculate and plot Ramachandran angles from a collection of
    structures
    """

    def __init__(
            self,
            input_structures,
            output_dir,
            plot_type,
            filter_residue_type=None,
            existing_rama_ditribution=None
        ) -> None:

        self.input_structures = input_structures
        self.output_dir = output_dir
        self.plot_type = plot_type


    def read_structure_files(self):
        """
        Read in PDB/mmCIF files from a directory
        """
        pass

    def filter_structure_residues(self):
        """
        OPTIONAL
        Filter structures based on the given plot-type
        """
        pass

    def calculate_dihedral_angles(self):
        """
        Calculate the phi and psi dihedral angles of the input structures
        """
        pass

    def save_dihedral_angles(self):
        """
        OPTIONAL
        Save the calculated dihedral angles to a CSV file
        """
        pass

    def plot_ramachandran(self, image_format=["png"]):
        """
        Plot the Ramachandran plot of the calculated dihedral angles
        """
        pass

    def open_background(self):
        """
        OPTIONAL
        Open the pre-rendered background image of the Ramachandran plot
        """
        pass

    def render_background(self, colour="Blues"):
        """
        OPTIONAL
        Render the background of the Ramachandran plot
        """
        pass











# Main function
def main(
        pdb,
        itmod,
        model_num,
        itchain,
        chain_num,
        plot_type,
        out_dir,
        save,
        file_type
    ):

    # Selecting user's desired Ramachandran plot
    options = ["All", "General", "Glycine", "Proline", "Pre-proline", "Ile-Val"]

    logger.info(f"Importing structure data {pdb}")

    userpdb_df = extract_dihedrals(
        pdb_file_name=pdb,
        chain_id=chain_num
    )
    # Remove invalid dihedral angles/angles from ligands or non-canonical residues
    userpdb_df = userpdb_df.dropna()


    # User input determines background
    plot_type = options[int(plot_type)]
    # Out file name
    plot_name = str(out_dir + '/' + pdb[:-4] + '_' + plot_type + "RamachandranPlot_tmp")



    ########################################################
    #                IMPORTING REFERENCE DATA               #

    logger.info("Importing Top8000 library")

    # Top8000 peptide dataset. Pre-analysed
    top8000_df = select_angles(
			pd.read_csv(
				"src/top8000_dihedrals.csv.gz",
				compression="gzip"
			),
            plot_type
        )



    ########################################################
    #                SELECTING RESIDUE TYPE DATA               #

    if plot_type == "All":
        pass

    else:
        userpdb_df = userpdb_df.loc[userpdb_df["type"] == plot_type]



    ########################################################
    #                SAVING USER DATA (optional)               #

    if save:
        csv_file_name = str(plot_name[:-4] + ".csv")
        logger.info(f"Saving CSV as {csv_file_name}")
        userpdb_df.to_csv(csv_file_name, index=False)


    logger.info("Calculating Ramachandran angles")


    ########################################################
    #                    PLOTTING DATA                       #

    # Open config file
    with open("src/plot_design_params.json") as f:
        plot_config = json.load(f)

    # Plotting background
    logger.info("Generating background of favoured regions")

	# Genertating background: region of favoured dihedral angles
    render_background(top8000_df, plot_type, plot_name, plot_config["background_colour"])

    # Plotting user's PDB dihedral angles
    logger.info("Plotting Ramachandran diagram")

    plt.style.use("seaborn-v0_8-poster")

    fig, ax = plt.subplots(1,1, figsize=plot_config["figure_size"], tight_layout=True)        # Defining plot area.

    # ADDING COUNTOURS - Comment this section out to remove contour lines from plot area.
    add_contour(
        ax,
        top8000_df,
        contour_level=plot_config["contour_level_inner"],
        line_colour=plot_config["contour_line_color_inner"]
    )
    add_contour(
        ax,
        top8000_df,
        contour_level=plot_config["contour_level_outer"],
        line_colour=plot_config["contour_line_color_outer"],
        contour_alpha=0.3
    )

    # ADDING FAVOURED RAMACHANDRAN REGION IMAGE TO BACKGROUND
    ax.imshow(
        plt.imread(str(plot_name + ".png")),
        extent=[-195, 195, -195, 195], zorder=1
    )

    # ADDING GRIDLINES
    add_grid_lines(ax)

    # PLOTTING USER'S DIHEDRAL ANGLE DATA
    ax.scatter(
        userpdb_df["phi"],
        userpdb_df["psi"],
        s=15,
        color=plot_config["data_point_colour"],
        zorder=4,
        linewidths=0.5,
        edgecolor=plot_config["data_point_edge_colour"]
    )

    # AXES AESTHETICS/FEATURES
    extra_formatting(ax)

    # SAVING RAMACHANDRAN PLOT AS PNG IMAGE
    logger.info("Saving plot")

    if file_type == "png":

        # ... as PNG
        plt.savefig(
            str(plot_name[:-4] + ".png"),
            dpi=plot_config["out_resolution"],
            bbox_inches=0,
            pad_inches=None
        )

    else:
        # ... as PDF
        plt.savefig(
            f"{plot_name[:-4]}.{file_type}",
            bbox_inches=0,
            pad_inches=None
        )

    os.remove(str(plot_name + '.png'))

    plt.close()

    logger.info(f"Done - Plot saved to{plot_name[:-4]}.{file_type}")
