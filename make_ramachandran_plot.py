import logging
import argparse
import pathlib

from src.ramachandran_plotter import RamachandranPlotter, render_background
from src.argument_parser import parse_user_args

import logging

if __name__ == "__main__":

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

    user_args = parse_user_args()

    # Parse structure files
    structures = {}
    for structure_file in user_args.input_file:
        key = structure_file[0]
        if len(structure_file) == 1:
            structures[key] = ['*']
        else:
            structures[key] = structure_file[1:]

    # Initialise
    my_rama_plotter = RamachandranPlotter(
        input_structures=structures,
        output_dir=pathlib.Path(user_args.out_dir),
        plot_type="all",
        # existing_rama_ditribution=None,
        # remove_outliers=False,
        # strictly_canonical=False,
        # exclude_additional_residues=[],
        # exclude_pre_proline=False,
        # config_file=None
    )

    # Extract structures from input files in structures dictionary
    my_rama_plotter.extract_structures()

    # Filter residues
    my_rama_plotter.filter_residues(
        # exclude_physicochem=[] # Could be your own data
        # exclude_additional_residues=[], # Could be your own data
    )

    # Calculate Ramachandran angles
    my_rama_plotter.calculate_ramachandran_angles()

    # Render background
    top8000_angles = my_rama_plotter.load_representative_angles(
        # Could be your own data
    )
    render_background(top8000_angles, "tmp/backgrounds/top8000.png")

    my_rama_plotter.plot_ramachandran(path_background="tmp/backgrounds/top8000.png")
