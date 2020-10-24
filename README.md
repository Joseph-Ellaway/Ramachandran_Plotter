# Ramachandran Plotter


## Requirements:

Python 3.8

Modules:
- Numpy
- Pandas
- Matplotlib
- Scipy
- cv2
- sys
- subprocess

Phenix Suite (install [here](https://www.phenix-online.org/download/))


## Example Ouput Figures

<!-- ![](http://example.com/image.jpg)
![](./image.jpg) -->
<!-- The above is commented out in HTML --> 

## Run Instructions

	python3 DihedralCalculator.py <file-name>
	python3 RamachandranPlot.py <file-name> <plot-code>

where file-name is your PDB file name (include path if necessary). 
	  plot-code can be any of the following integers to determine the type of output plot desired:

	0 = All angles 
 	1 = General angles	(All angles excluding Gly, Pro, Ile and Val)
 	2 = Glycine angles 
 	3 = Trans-proline angles only 
 	4 = Cis-proline angles only 
 	5 = Proline (cis/trans) angles only 
 	6 = Pre-proline angles only (angles preceding a Pro residue)
 	7 = Ile-Val angles only 

The first command generates a CSV of your input peptides dihedral angles using the phenix.ramalyze function included in the Phenix Suite of 
structural tools. 

The second command takes the CSV (where multiple peptides have had their dihedral angles calculated, the peptide PDB file name must also be given)
and creates a Ramachandran plot. The type of plot returned is determined by the plot code. 

Backgrounds to Ramachandran plots are generated using the results of phenix.ramalyze applied to all peptides from the Top8000 peptide database. 
These are peptides for which models have been solved at very high resolutions and dihedral angles are assumed to be at their true values. 

Several parameters can be easily adjusted to change the appearence of the returned graph. 