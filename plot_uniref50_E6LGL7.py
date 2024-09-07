import logging
import pathlib

from src.ramachandran_plotter import RamachandranPlotter, render_background

import logging

if __name__ == "__main__":

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

    # Input structures
    base_path = pathlib.Path("example_data", "uniref50_E6LGL7")
    structures = {
        base_path / "5dvy_updated.cif" : ['*'],
        base_path / "6g88_updated.cif" : ['*'],
        base_path / "6mki_updated.cif" : ['*'],
        base_path / "8f3i_updated.cif" : ['*'],
        base_path / "8f3o_updated.cif" : ['*'],
        base_path / "8f3t_updated.cif" : ['*'],
        base_path / "8f3y_updated.cif" : ['*'],
        base_path / "5e31_updated.cif" : ['*'],
        base_path / "6mka_updated.cif" : ['*'],
        base_path / "6mkj_updated.cif" : ['*'],
        base_path / "8f3j_updated.cif" : ['*'],
        base_path / "8f3p_updated.cif" : ['*'],
        base_path / "8f3u_updated.cif" : ['*'],
        base_path / "8f3z_updated.cif" : ['*'],
        base_path / "6bsq_updated.cif" : ['*'],
        base_path / "6mkf_updated.cif" : ['*'],
        base_path / "8f3f_updated.cif" : ['*'],
        base_path / "8f3l_updated.cif" : ['*'],
        base_path / "8f3q_updated.cif" : ['*'],
        base_path / "8f3v_updated.cif" : ['*'],
        base_path / "8f67_updated.cif" : ['*'],
        base_path / "6bsr_updated.cif" : ['*'],
        base_path / "6mkg_updated.cif" : ['*'],
        base_path / "8f3g_updated.cif" : ['*'],
        base_path / "8f3m_updated.cif" : ['*'],
        base_path / "8f3r_updated.cif" : ['*'],
        base_path / "8f3w_updated.cif" : ['*'],
        base_path / "8u55_updated.cif" : ['*'],
        base_path / "6g0k_updated.cif" : ['*'],
        base_path / "6mkh_updated.cif" : ['*'],
        base_path / "8f3h_updated.cif" : ['*'],
        base_path / "8f3n_updated.cif" : ['*'],
        base_path / "8f3s_updated.cif" : ['*'],
        base_path / "8f3x_updated.cif" : ['*'],
    }

    # Initialise
    my_rama_plotter = RamachandranPlotter(
        input_structures=structures,
        output_dir=pathlib.Path("example_data", "example_output"),
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
