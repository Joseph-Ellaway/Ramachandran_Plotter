# Ramachandran Plotter (v3.0.0)


## Requirements:

Python 3.11

Modules:
- Numpy
- Pandas
- Matplotlib
- Scipy
- CV2
- Biopython
- OS
- Argparse

Install with pip:

	pip install -r requirements


## New config file

// Output figure size in inches.
// Percentile of dihedral angles for inner contour lines (e.g. contour_level=96 means the area bounded by the contour line represents the range of angles in which 96% of all dihedral from the Top800 peptide DB fall within).
// Percentile of dihedral angles for outer contour lines
// Inner contour line colour.
// Colour of outer contour lines
// Output figure resolution. Not required if saving file as PDF
// Colour of data points for each Phi-Psi dihedral angle pair.
// Colour of data point's border.
// Colour map of background plot. Refer to https://matplotlib.org/stable/tutorials/colors/colormaps.html for colormap options

## Run Instructions: CLI

The script `make_ramachandran_plot.py` may be used to access the `RamachandranPlotter()`
object, without needing to write any code. This option will be preferred by users who
need a plot, given a set of structures. There are several optional arguments, and
therefore different ways of running the `make_ramachandran_plot.py` script. The
following are examples in order of increasing complexity.

### 1) Minimal run example

The simplest way to run the script is as follows


```shell
python3 make_ramachandran_plot.py \
	--input_file /path_to_file/<file-name-1.cif> A C \
	--input_file /path_to_file/<file-name-2.cif> A B \
	--input_file ...
	--input_file /path_to_file/<file-name-n.cif> * \

```

The input to `--input_file` is formatted as the path to your mmCIF (or PDB) file,
followed by a space-separated list of chain IDs. By default, author-specified chain IDs
(`auth_asym_id`) are used, but this can be overridden to use `struct_asym_id` by parsing
in the additional argument `--use-struct-asym-ids` or `-u`. To parse all chains in the
file, use the wildcard `*`.

For example:

```shell
python3 make_ramachandran_plot.py --input_file example_data/uniref50_E6LGL7/5e31_updated.cif *
```



Optional arguments:

	--help				: Prints summary of arguments
	--verbose			: Increase output verbosity
	--models <int>		: Desired model number (default = use all models). Model number corresponds to order in PDB file.
	--chains <int>		: Desired chain number (default = use all chains). Chain number corresponds to order in PDB file.
	--out_dir <path>	: Out directory. Must be available before-hand.
	--plot_type <int>	: Type of angles plotted Ramachandran diagram. Options detailed below.
	--save_csv			: Saves calculated dihedral angles in a separate CSV file.

```--plot_type <int>``` can be any of the following integers to determine the type of output plot desired:

	0 	: All angles
	1 	: General (All residues bar Gly, Pro, Ile, Val and pre-Pro)
	2 	: Glycine
	3 	: Proline (cis and trans)
	4 	: Pre-proline (residues preceeding a proline)
	5 	: Ile or Val

Backgrounds to Ramachandran plots are generated using dihedral angle data from peptide structures solved at high resolution from the Top8000 peptide database.

These are peptides for which models have been solved at very high resolutions and dihedral angles are assumed to be at their true values.

Several parameters can be easily adjusted to change the appearance of the returned graph.

#### All angle plot

	python RamachandranPlotter.py --pdb 6GVE.pdb --plot_type 0

<img src="./example_plot/AllRamachandranPlot.png" style="zoom:35%;" />

## Adjustable Variables (recommended)

| Parameter | Variable name | Description |
| :--- | :--- | :--- |
|Figure size| **```figure_size```** |Adjusts the output figure size (inches) as a tuple|
|Figure resolution | **```out_resolution```**|Output final figure resolution (high values will slow the process down)|
|Inner contour line level |**```contour_level_inner```**|Level at which to draw the inner contour lines. Should be a value between 0-100 to represent to percemtile at which dihedral angles from the Top8000 peptide DB fall within. e.g. a value 96 coresponds to the area where at least 96 % of dihedral angles fall within. |
| Outer contour line level | **```contour_level_outer```** |Level at which to draw the outer contour lines. Should be a value between 0-100 to represent to percemtile at which dihedral angles from the Top8000 peptide DB fall within. e.g. a value 15 coresponds to the area where at least 15 % of dihedral angles fall within. |
|Favoured region colour |**```background_colour```**| Colour of the favoured dihedral angle region data points are plotted against*. |
|Inner contour line colour |**```contour_line_color_inner```**|Colour of inner contour lines. |
|Outer contour line colour |**```contour_line_color_outer```**|Colour of outer contour lines. |
|Data point colour |**```data_point_colour```** |Colour of data point for all dihedral angle pairs. |
|Data point edge colour |**```data_point_edge_colour```** |Colour of the borders for data point colours.|

\* options for sequential colour maps (recommended):
[````'Greys'````, ````'Purples'````, ````'Blues'````, ````'Greens'````, ````'Oranges'````, ````'Reds'````,
````'YlOrBr'````, ````'YlOrRd'````, ````'OrRd'````, ````'PuRd'````, ````'RdPu'````, ````'BuPu'````,
````'GnBu'````, ````'PuBu'````, ````'YlGnBu'````, ````'PuBuGn'````, ````'BuGn'````, ````'YlGn'````]

Find a complete description of available colour maps and how to make your own at: <a href="https://matplotlib.org/stable/tutorials/colors/colormaps.html">https://matplotlib.org/stable/tutorials/colors/colormaps.html</a>

## Updates for version 2.0.2
Only one command is required to generate Ramachandran plot. Before, the CSV of dihedral angles had to be generated using a separate script. In v2.0.2, a CSV of dihedral angles can be generated by parsing the ```--save_csv``` argument.

The dependence previous versions of Ramachandran-Plotter had on Phenix has been replaced with Biopython.

**01/08/2021**: User can specify desired file type of final plot on execution of ```RamachandranPlotter.py```. Options: ```PNG``` (default, 96 dpi. Refer to table above to change), ```PDF```, ```SVG```, ```EPS``` or ```PS```.



## Bug fixes:

### 07/10/2023

Crash fixes:
 - Dataframe.append(...) is deprecated since pandas 2.0, replaced by pandas.concat(...)
 - Matplotlib style *seaborn-poster* was renamed *seaborn-v0_8-poster* since 3.6.3 version

Bug fixes:
 - Temporary png image was not deleted on Windows, replace *os.command()* by *os.remove()*

### 27/03/2022

- Improved readability:
  - Better comment brevity
  - Majority of code fit within column width of 88, where possible/appropriate

Crash fixes:

- Replaced LaTeX package implementation from Matplotlib (for rendering phi and psi symbols) with unicode chars.

### 01/08/2021

Crash fixes:
 - Top8000 file decompression error has been fixed.
 - Plotting Proline-only plots has been patched.
 - Corrections to dependency install instructions.

Bug fixes:
 - Contour lines for specified plot now render correctly.
 - Background plot of all dihedral angles is now only created with ```plot_type 0```, if specified.
 - Out directory bug fixed.
 - PDB file name is now appended to out files.

 Optimisation improvements:
 - Specified angles, if user does not want to plot all dihedral types, are now selected once at the start.



### Author Details

Joseph I. J. Ellaway

josephellaway@gmail.com

MSc Bioinformatics and Theoretical Systems Biology, Imperial College London
