import matplotlib.pyplot as plt
import pandas as pd
import logging
import json
import pathlib
import gemmi
import tempfile

# Package functions
from .rama_calc_utils import *
from .plot_utils import *
from .argument_parser import *

# Logging
logger = logging.getLogger(__name__)

def render_background(
        background_angles: np.array,
        path_save: pathlib.Path|str = None,
        **kwargs
    ):
    """
    TODO: this needs cleaning up
    Render the background of the Ramachandran plot
    """
    logger.info("Rendering background")

    # Optional parameters
    background_colour = kwargs.get("background_colour", "Blues")

    # Plot angles to histogram
    render_fig_size = (10,10)
    background_fig, background_ax = plt.subplots(
        1,1,
        figsize=render_fig_size,
        tight_layout=True
    )
    background_fig.patch.set_visible(False)
    background_ax.set_axis_off()
    hist_kwargs = {
        "bins" : 140,
        "cmap" : background_colour,
        "norm" : colors.PowerNorm(0.1),
        "alpha" : 1
    }

    background_ax.hist2d(
        background_angles[:,0],
        background_angles[:,1],
        **hist_kwargs
    )

    tmp_file_name = f"{tempfile.NamedTemporaryFile().name}.png"
    logger.debug(f"Temporary file name: {tmp_file_name}")
    save_and_close_fig(tmp_file_name, 80)

    # Apply post-processing to the background image
    background_fig, background_ax = plt.subplots(
        1,1,
        figsize=render_fig_size,
        tight_layout=True
    )
    background_fig.patch.set_visible(False)
    img_wkargs = {
        "cmap" : str(background_colour + "_r"),
        "alpha" : 1
    }
    background_ax.set_axis_off()
    rama_plot = cv2.imread(tmp_file_name, 0)
    blurred_rama_plot = ndimage.gaussian_filter(rama_plot, sigma=0.3)
    smoothed_rama_plot = ndimage.percentile_filter(
        blurred_rama_plot,
        percentile=90,
        size=20
    )
    background_ax.imshow(smoothed_rama_plot, **img_wkargs)

    save_and_close_fig(path_save, 96)
    logger.info("Background rendered")


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


def save_and_close_fig(out_file_name, resolution):			# dtypes : string, int
	"""
	Saves the plot area (defined before function called) as a PNG image to given file
	name and resolution.
	"""
	plt.savefig(out_file_name, dpi=resolution, bbox_inches=0, pad_inches=None)
	plt.close()




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
        self.output_dir = pathlib.Path(pathlib.Path(output_dir).stem)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.output_plot_fname = pathlib.Path(output_dir)
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


    def extract_structures(self):
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


    def filter_residues(
            self,
            exclude_physicochem: list = None,
            exclude_additional_residues: list = None
        ):
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

        if exclude_physicochem:
            for property in exclude_physicochem:
                match property:
                    case "hydrophobic":
                        self.excluded_residues.extend(
                            ["ILE", "VAL", "LEU", "MET", "PHE", "TRP"]
                        )
                    case "polar":
                        self.excluded_residues.extend(
                            ["SER", "THR", "ASN", "GLN", "CYS", "TYR"]
                        )
                    case "charged":
                        self.excluded_residues.extend(
                            ["ARG", "LYS", "ASP", "GLU", "HIS"]
                        )
                    case "aromatic":
                        self.excluded_residues.extend(
                            ["PHE", "TRP", "TYR"]
                        )
                    case "positive":
                        self.excluded_residues.extend(
                            ["ARG", "LYS", "HIS"]
                        )
                    case "negative":
                        self.excluded_residues.extend(
                            ["ASP", "GLU"]
                        )
                    case "small":
                        self.excluded_residues.extend(
                            ["ALA", "GLY", "SER", "THR", "CYS"]
                        )
                    case "tiny":
                        self.excluded_residues.extend(
                            ["ALA", "GLY", "SER"]
                        )
                    case "large":
                        self.excluded_residues.extend(
                            ["ILE", "LEU", "MET", "PHE", "TRP", "TYR"]
                        )
                    case _:
                        logger.info(f"Unknown physico-chemical property: {property}")

        # Optional additional residues to exclude
        if exclude_additional_residues:
            self.excluded_residues.extend(exclude_additional_residues)

        # Convert to set for faster lookup
        self.excluded_residues = set(self.excluded_residues)



    def calculate_ramachandran_angles(self):
        """
        Calculate the phi and psi dihedral angles of the input structures
        """

        # List to store the calculated angles
        self.all_input_angles = [
            # [phi, psi] -- 2D array
        ]
        # angles = {}

        for label, chain in self.polymers.items():
            logger.info(f"Calculating Ramachandran angles for {label}")
            for res in chain.get_polymer():
                if res.name in self.excluded_residues:
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
                self.all_input_angles.append(angle)

        # Nothing calculated
        if len(self.all_input_angles) == 0:
            logger.error("No angles calculated")
            return None

        # Convert to numpy array and degrees
        self.all_input_angles = np.array(self.all_input_angles)
        self.all_input_angles = np.degrees(self.all_input_angles)

        # Optional
        return self.all_input_angles


    def load_representative_angles(
            self,
            path_to_data: pathlib.Path|str = "src/top8000_dihedrals.csv.gz"
        ):
        """
        OPTIONAL
        Load the pre-calculated representative Ramachandran angles
        TODO: Consider putting some filters in place
        """
        logger.info("Loading representative Ramachandran angles")

        repr_table = pd.read_csv(
            path_to_data,
            compression="gzip"
        )

        # Select only phi and psi angles not in excluded residues
        self.repr_angles = repr_table[["phi", "psi"]][
            ~repr_table["residueID"].isin(self.excluded_residues)
        ].to_numpy()

        return self.repr_angles






    def plot_ramachandran(
            self,
            image_formats=["png",],
            path_background=None
        ):
        """
        Plot the Ramachandran plot of the calculated dihedral angles
        """

        # Render the background image
        if not path_background:
            path_background = f"tmp/backgrounds/background.png"
            render_background(
                self.repr_angles,
                path_save=path_background
            )


        plot_config = json.load(open("src/plot_design_params.json"))

        config_file = pathlib.Path("plot_design_params.json")
        fig, ax = plt.subplots(1,1, figsize=plot_config["figure_size"], tight_layout=True)        # Defining plot area.

        # Add contour lines
        add_contour(
            ax,
            self.repr_angles,
            contour_level=plot_config["contour_level_inner"],
            line_colour=plot_config["contour_line_color_inner"]
        )
        add_contour(
            ax,
            self.repr_angles,
            contour_level=plot_config["contour_level_outer"],
            line_colour=plot_config["contour_line_color_outer"],
            contour_alpha=0.3
        )

        # ADDING FAVOURED RAMACHANDRAN REGION IMAGE TO BACKGROUND
        # TODO: Update this to dynamic name
        ax.imshow(
            plt.imread(path_background),
            extent=[-195, 195, -195, 195], zorder=1
        )

        # ADDING GRIDLINES
        add_grid_lines(ax)

        # PLOTTING USER'S DIHEDRAL ANGLE DATA
        ax.scatter(
            self.all_input_angles[:,0],     # phi
            self.all_input_angles[:,1],     # psi
            s=15,
            color=plot_config["data_point_colour"],
            zorder=4,
            linewidths=0.5,
            edgecolor=plot_config["data_point_edge_colour"]
        )

        # AXES AESTHETICS/FEATURES
        limits = (-180, 180)

        ax.set_xlim(limits)
        ax.set_ylim(limits)
        ax.set_xlabel(u"\u03A6 (\u00B0)")	# phi
        ax.set_ylabel(u"\u03A8 (\u00B0)")	# psi
        ax_linewidth = 2
        ax.spines["left"].set_linewidth(ax_linewidth)
        ax.spines["bottom"].set_linewidth(ax_linewidth)

        # SAVING RAMACHANDRAN PLOT AS PNG IMAGE
        logger.info("Saving plot")
        for image_format in image_formats:
            plot_dir_fname = self.output_plot_fname
            if image_format == "png":
                # ... as PNG
                plt.savefig(
                    plot_dir_fname,
                    dpi=plot_config["out_resolution"],
                    bbox_inches=0,
                    pad_inches=None
                )
            else:
                plt.savefig(
                    plot_dir_fname,
                    bbox_inches=0,
                    pad_inches=None
                )
            logger.info(f"Plot saved as {plot_dir_fname}")

            plt.close()



    def save_ramachandran_angles(self, fname: str = "ramachandran_angles.csv"):
        """
        Save the calculated dihedral angles to a CSV file
        """
        np.savez_compressed(fname, angles=self.all_input_angles)
        logger.info(f"Ramachandran angles saved to {fname}")
