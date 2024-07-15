import os
import matplotlib.pyplot as plt
import pandas as pd
import logging
import json
import pathlib
import gemmi
import sys
from math import degrees
import gemmi

# Package functions
from .rama_calc_utils import *
from .plot_utils import *
from .argument_parser import *

# Logging
logger = logging.getLogger(__name__)



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




def calculate_rama_angles(pdb_file, chain_id=None):
    """
    Calculate the phi and psi dihedral angles of a PDB file
    """
    # # Extract dihedral angles
    # dihedral_angles = extract_dihedrals(
    #     pdb_file_name=pdb_file,
    #     chain_id=chain_id
    # )

    # # Remove invalid dihedral angles/angles from ligands or non-canonical residues
    # dihedral_angles = dihedral_angles.dropna()

    # return dihedral_angles



























def map_auth_to_label_asym_ids(chain):
    """
    Convert auth_seq_id to label_seq_id
    """
    # Load the structure
    structure = gemmi.read_structure("6gve_updated.cif")

    # Get the first model
    model = structure[0]
    print(dir(model))


    for chain in model:
        # print(dir(chain))
        print(
            chain.name,
            chain.subchains()[0].subchain_id(),
        )


























class RamachandranPlotter:
    """
    Object to parse, calculate and plot Ramachandran angles from a collection of
    structures
    """

    def __init__(
            self,
            input_structures: dict,
            output_dir : pathlib.Path|str,
            plot_type : str,
            filter_residue_type=None,
            existing_rama_ditribution=None,
            remove_outliers=False,
            strictly_canonical=False,
            exclude_additional_residues=[],
            exclude_pre_proline=False,
            config_file=None        # Specify default
        ) -> None:
        """
        Example formats:

        input_structures = {
            "file1.pdb" : ["A", "B", "C"],    # PDB format, chain IDs specified
            "file2.cif" : ["A", "B", "C"],    # mmCIF format, chain IDs specified
            "file3.pdb" : ['*'],              # PDB format, use all chains
            "file4.cif" : ['*'],              # mmCIF format, use all chains
        }
        """
        self.input_structures = input_structures
        self.output_dir = output_dir
        self.plot_type = plot_type

        self.excluded_residues = exclude_additional_residues
        self.exclude_pre_proline = exclude_pre_proline
        self.strictly_canonical = strictly_canonical
        self.canonical_residues = (
            "LEU", "ALA", "GLY", "VAL", "GLU", "SER", "LYS", "ASP", "THR", "ILE",
            "ARG", "PRO", "ASN", "PHE", "GLN", "TYR", "HIS", "MET", "CYS", "TRP"
        )

    def convert_auth_to_label_asym_ids(self):
        """
        OPTIONAL
        Convert auth_seq_id to label_seq_id
        """
        pass


    def read_structure_files(self):
        """
        Read in PDB/mmCIF files from a directory
        """

        self.polymers = {}

        # NOTE: Assume auth_asym_id for now
        for structure_path, chains in self.input_structures.items():
            # Read in the structure
            structure = gemmi.read_structure(str(structure_path))
            model = structure[0]

            # All
            if chains == ['*']:
                chains = [chain.name for chain in model]

            for chain in chains:
                label = f"{structure.name}_{chain}"
                residue_span = model[chain]
                self.polymers[label] = residue_span

        # Optional
        return self.polymers


    def filter_structure_residues(self, exclude_additional_residues: list = None):
        """
        OPTIONAL
        Filter structures based on the given plot-type
        """

        # Defaults
        match self.plot_type:
            case "all":
                pass
            case "general":
                self.excluded_residues = ["GLY", "PRO", "ILE", "VAL"]
            case "glycine":
                self.excluded_residues = ["GLY"]
            case "proline":
                self.excluded_residues = ["PRO"]
            case "pre-proline":
                self.exclude_pre_proline = True
            case "ile-val":
                self.excluded_residues = ["ILE", "VAL"]
            case _:
                logger.info("No plot type specified. Defaulting to all residues")

        # Optional additional residues to exclude
        if exclude_additional_residues:
            self.excluded_residues.extend(exclude_additional_residues)



    def calculate_ramachandran_angles(self):
        """
        Calculate the phi and psi dihedral angles of the input structures
        """

        # List to store the calculated angles
        angles = []

        # structure = gemmi.read_structure("6gve_updated.cif")

        for label, chain in self.polymers.items():
            logger.info(f"Calculating Ramachandran angles for {label}")
            for res in chain.get_polymer():
                if res in self.excluded_residues:
                    # Exluded residue
                    continue

                prev_res = chain.previous_residue(res)
                next_res = chain.next_residue(res)

                if not prev_res and not next_res:
                    # No previous or next residue
                    logger.debug(f"Residue {res} has no previous or next residue")
                    continue

                if self.exclude_pre_proline and prev_res.name == 'PRO':
                    # Exclude pre-proline residues (conditional)
                    logger.debug(f"Excluding pre-proline residue for {res}")
                    continue

                if self.strictly_canonical and res.name not in self.canonical_residues:
                    # Only canonical residues
                    logger.debug(
                        f"Excluding non-canonical residue {res}, name: {res.name}"
                    )
                    continue

                # Calculate phi and psi angles
                angle = gemmi.calculate_phi_psi(prev_res, res, next_res)
                angles.append(angle)

        angles = np.array(angles)
        angles = np.degrees(angles)

        print(angles)

        # Optional
        return angles



    def save_dihedral_angles(self):
        """
        OPTIONAL
        Save the calculated dihedral angles to a CSV file
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

    def plot_ramachandran(
            self,
            image_format=["png"],
            save=True,
            config_file=pathlib.Path("plot_design_params.json")
        ):
        """
        Plot the Ramachandran plot of the calculated dihedral angles
        """
        pass
