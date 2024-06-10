


import logging
import argparse

from src.ramachandran_plotter import main
from src.argument_parser import collect_user_args


if __name__ == "__main__":

    # Loading user's input arguments
    (
        pdb,
        itmod,
        model_num,
        itchain,
        chain_num,
        plot_type,
        out_dir,
        verb,
        save,
        file_type
    ) = collect_user_args()


    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )
    logging.debug("Starting Ramachandran Plotter")

    main(
        pdb,
        itmod,
        model_num,
        itchain,
        chain_num,
        plot_type,
        out_dir,
        save,
        file_type
    )
