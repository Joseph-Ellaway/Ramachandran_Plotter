
''' ================================================================================================================

	Functions here are called by main() to created a Pandas DataFrame of dihedral angle information from a given PDB
	file. Functions can also be used elsewhere for standalone use. 
	
	Version 2.0.1:
	 - Relies on the easily accessible Biopython package, rather than Phenix as in versions <2.0 (as facilitated by 
	 functions here)
	 - User arguments can be now easily parsed in from the command line
	 - If required, the script could be implemented into existing protein analysis pipelines by importing this 
	 function ( main() )

	Author information:
	 - Joseph I. J. Ellaway
	 - josephellaway@gmail.com
	 - https://github.com/Joseph-Ellaway

	================================================================================================================ '''

# Removes warning messages from Biopython. Comment out to debug
import warnings
warnings.filterwarnings("ignore")

import Bio.PDB
import pandas as pd
import numpy as np
import math


import warnings
warnings.filterwarnings("ignore")

import Bio.PDB

import pandas as pd
import numpy as np
import math


<<<<<<< HEAD
def ResidueNames(chain):
	'''
	===========================================================================================================
	Takes a Biopython PDB polypeptide object and returns a two lists: residue names and their position indices
	===========================================================================================================
	'''

	# Initialising output lists
	res_names = []
	res_indices = []

	# Iterate over residues in chain
	for index, residue in enumerate(chain):

		res_names.append(residue.resname) 	# Residue name (3-letter code)
		res_indices.append(residue.id[1])	# Residue position in chain

	return res_names, res_indices



def ToDegrees(radian_list):
	'''
	===============================================
	Converts a list of floats in radians to degrees
	===============================================
	'''

	# Initialise output list
	out_list = []

	# Iterate over floats (in radians in list) and convert each to degrees
	for i in radian_list:
		try:
			degrees = math.degrees(i)
			out_list.append(degrees)
		except:
			out_list.append(np.nan)  # Handles NaN values

	return out_list



def CalcDihedrals(polypep):
	'''
	=========================================================================================
	Takes a Biopython polypeptide object and returns two lists of: Phi angles and Psi angles
	=========================================================================================
	'''

	# Calculate dihedral angles for residues in polypeptide
	angles = Bio.PDB.Polypeptide.Polypeptide(polypep)
	angles = angles.get_phi_psi_list()

	# Assign phi and psi angles to separate list variables
	phis = []
	psis = []

	for angle_pair in angles:
		phis.append(angle_pair[0])
		psis.append(angle_pair[1])

	phis = ToDegrees(phis)
	psis = ToDegrees(psis)

	# Return phi and psi angles as separate list variables
	return phis, psis



def AminoAcidType(residue_names):
	'''
	=========================================================================================================
	Takes a list of residue names (3-letter codes) and classifies them into one of six categories:
		- Glycine
		- Proline
		- Isoleucine/valine
		- Pre-proline (any residue preceeding a Pro)
		- General (any canonical residue that is not classified by the above)
	Outputs a list of class types, length is equal to input list. Invalid residues are classed as NaN. 
	=========================================================================================================
	'''

	# Initialise output list
	residue_types = []

	count = 0 		# Iterable

	# Iterate over residue names in input list
	while count < len(residue_names):

		aa = residue_names[count]

		# Sometimes residue codes include prefix (for details on their PTM). They must be removed here.
		if len(aa) > 3:
			aa = aa[:-3].upper()
		else:
			aa = aa.upper()

		# Check if the residue is a canonical amino acid -- entries could be ligands and should not be processed. 
		if aa in ['MET', 'SER', 'ASN', 'LEU', 'GLU', 'LYS', 'GLN', 'ILE', 'ALA', 'ARG', 'HIS', 'CYS', 'ASP', 'THR', 'GLY', 'TRP', 'PHE', 'TYR', 'PRO', 'VAL']:

			# Glycine check
			if aa == 'GLY':
				aa_type = 'Glycine'

			elif aa == 'PRO':
				aa_type = "Proline"

			# Isoleucine/valine check
			elif aa == 'ILE' or aa == 'VAL':
				aa_type = 'Ile-Val'

			# General check
			elif aa not in ['GLY', 'PRO', 'ILE', 'VAL']:
				aa_type = 'General'
			
			# Entry is invalid 
			else:
				aa_type = np.nan

			# Pre-proline check
			try:
				if residue_names[count + 1] == 'PRO':
					aa_type = 'Pre-proline'
				else:
					pass
			except:
				pass
		else:
			aa_type = np.nan
			

		residue_types.append(aa_type)
		count += 1

	return residue_types



def ChainSummary(polypep):
	'''
	===========================================================================================
	Returns relevant information on a Biopython polypeptide object for downstream processing:
		- Phi/Psi angles
		- Residue names/position indices
		- Residue type
		- Chain ID
	... in a Pandas DataFrame.
	===========================================================================================
	'''

	# Calculate dihedral angles for residues in chain and add them to separate list variables
	chain_phis, chain_psis = CalcDihedrals(polypep)

	# Return residue names and position indices within polypeptide chain
	chain_resnames, chain_resindices = ResidueNames(polypep) 

	# Return the type of the amino acid
	chain_types = AminoAcidType(chain_resnames)
=======
########################################################
#						FUNCTIONS 					   #
########################################################



def ResidueNames(chain):

	res_names = []
	res_indices = []

	for index, residue in enumerate(chain):

		res_names.append(residue.resname)
		res_indices.append(residue.id[1])

	return res_names, res_indices



def ToDegrees(radian_list):

	out_list = []

	for i in radian_list:
		try:
			degrees = math.degrees(i)
			out_list.append(degrees)
		except:
			out_list.append(np.nan)

	return out_list




def CalcDihedrals(polypep):

	# Calculate dihedral angles for residues in polypeptide
	angles = Bio.PDB.Polypeptide.Polypeptide(polypep)
	angles = angles.get_phi_psi_list()

	# Assign phi and psi angles to separate list variables
	phis = []
	psis = []

	for angle_pair in angles:
		phis.append(angle_pair[0])
		psis.append(angle_pair[1])

	phis = ToDegrees(phis)
	psis = ToDegrees(psis)

	# Return phi and psi angles as separate list variables
	return phis, psis

def AminoAcidType(residue_names):

	residue_types = []

	count = 0

	while count < len(residue_names):

		aa = residue_names[count]

		# Sometimes residue names include prefix for details on their PTM 
		if len(aa) > 3:
			aa = aa[:-3]

		else:
			pass

		# Check if the residue is a canonical amino acid -- entries could be ligands and should not be processed. 
		if aa in ['MET', 'SER', 'ASN', 'LEU', 'GLU', 'LYS', 'GLN', 'ILE', 'ALA', 'ARG', 'HIS', 'CYS', 'ASP', 'THR', 'GLY', 'TRP', 'PHE', 'TYR', 'PRO', 'VAL']:

			try:
				if residue_names[count + 1] == 'Pro':
					aa_type = 'Pre-proline'
				else:
					pass
			except:
				pass

			if aa == 'GLY':
				aa_type = 'Glycine'

			elif aa == 'ILE' or aa == 'VAL':
				aa_type = 'Ile-Val'

			elif aa not in ['GLY', 'PRO', 'ILE', 'VAL']:
				aa_type = 'General'
			
			else:
				aa_type = np.nan
		else:
			aa_type = np.nan
			

		residue_types.append(aa_type)
		count += 1

	return residue_types

def ChainSummary(polypep):
	# Calculate dihedral angles for residues in chain and add them to separate list variables
	chain_phis, chain_psis = CalcDihedrals(polypep)     ## Done

	# Return residue names and position indices within polypeptide chain
	chain_resnames, chain_resindices = ResidueNames(polypep)   ## Done

	# Return the type of the amino acid
	chain_types = AminoAcidType(chain_resnames)

	# Generate chainID columns
	chain_IDs = [polypep.id] * len(chain_phis)

	# Add all polypep information to model DataFrame
	chain_summary = {
	'chainID' : chain_IDs,
	'residueName' : chain_resnames,
	'residueIndex' : chain_resindices,
	'phi' : chain_phis,
	'psi' : chain_psis,
	'type': chain_types
	}

	chain_summary = pd.DataFrame.from_dict(chain_summary)

	return chain_summary



def ModelDihedrals(model, model_num, iter_chains=True, chain_id=None):
	'''
	=======================================================
	Generates a Pandas DataFrame of phi/psi angles (and other information) from a given PDB model 
	=======================================================
	'''

	model_summaryDF = pd.DataFrame(columns=['chainID','residueName','residueIndex','phi','psi','type'])	

	# Iterate over all chains in model
	if iter_chains:

		for chain in model:
			chain_summaryDF = ChainSummary(chain)
			model_summaryDF = model_summaryDF.append(chain_summaryDF, ignore_index=True)

	# If multiple chains are present in the structure, default is to calculate dihedrals from all chains
	else:
		chain = model[chain_id]
		chain_summaryDF = ChainSummary(chain)
		model_summaryDF = model_summaryDF.append(chain_summaryDF, ignore_index=True)

	# Append model number information to final DataFrame
	model_ID_list = [model_num] * len(model_summaryDF)
	model_summaryDF.insert(loc=0, column='ModelID', value=model_ID_list)

	return model_summaryDF




def ExtractDihedrals(pdb_file_name=None, iter_models=True, model_number=0, iter_chains=True, chain_id=None):
	'''
	==============================================================================================================
	Generates a Pandas DataFrame of phi/psi angles (and other information) from a given PDB model 
	==============================================================================================================
	'''

	# Checks if a PDB file name is given
	if pdb_file_name:

		try:
			pdb_code = pdb_file_name[:-4]
			pdb_summaryDF = pd.DataFrame(columns=['ModelID','chainID','residueName','residueIndex','phi','psi','type'])

			if iter_models:

				models = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_file_name)

				for model in models:

					model_dihedrals = ModelDihedrals(model, model_number, iter_chains, chain_id)
					pdb_summaryDF = pdb_summaryDF.append(model_dihedrals, ignore_index=True)

					model_number += 1

			else:
				try:
					model = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_file_name)[model_number]
					model_dihedrals = ModelDihedrals(model, model_number)
					pdb_summaryDF = pdb_summaryDF.append(model_dihedrals, ignore_index=True)
				except:
					print('\n \n  !!! Invalid model number entered !!! \n \n')

			# Append PDB code information to final DataFrame
			pdb_list = [pdb_code] * len(pdb_summaryDF)
			pdb_summaryDF.insert(loc=0, column='PDBCode', value=pdb_list)

			return pdb_summaryDF

		except:
			print('\n \n  !! Invalid PDB file !! \n \n ' )

	else:
		print('\n \n !! No PDB file specified !! \n')
		print(' Specify PDB file using: \n \n     --PDB /path_to_file/<filename.pdb> \n \n')



# Testing
# output_DF = ExtractDihedrals('6GVE.pdb')

# aa_type_check = AminoAcidType(['MET', 'SER', 'ASN', 'LEU', 'GLU', 'LYS', 'GLN', 'ILE', 'ALA', 'ARG', 'HIS', 'CYS', 'ASP', 'THR', 'GLY', 'TRP', 'PHE', 'TYR', 'PRO', 'VAL'])
# print(aa_type_check)

# print(output_DF)







# def RamalyzeParseToCSV(stdout_list):
# 	'''
# 	=======================================================
# 	Takes every element from the penix.ramalyse results list and parse them into a .csv file to be 
# 	read with Pandas in future graph making scripts.
# 	=======================================================
# 	'''
# 	out_csv_name = str(sys.argv[1][:-4] + '_DihedralAngles.csv')
# 	ramalyze_results = open(out_csv_name, "w")
# 	ramalyze_results.write(str('chainID,residue,residueID,score%,phi,psi,evaluation,type' + '\n'))		# Adds the PDB file name to the start of the description row. 
# 	ramalyze_results_listolists = []

# 	for residue_info in stdout_list:
# 		residue_info_split = list(residue_info.split(' '))
# 		ramalyze_results_listolists.append(residue_info_split)

# 	ramalyze_results_listolists_no_spaces = []

# 	for embedded_list in ramalyze_results_listolists:
# 		rama_info_no_spaces = []
# 		for item in embedded_list:
# 			if item != '':
# 				rama_info_no_spaces.append(item)
# 			else:
# 				pass

# 		ramalyze_results_listolists_no_spaces.append(rama_info_no_spaces)

# 	for rama_info_list_no_spaces in ramalyze_results_listolists_no_spaces:
# 		chainID = rama_info_list_no_spaces[0]
# 		res_num = rama_info_list_no_spaces[1]
# 		the_good_stuff = max(rama_info_list_no_spaces, key=len)		# len key used so to avoid getting high values, such as 'v's, in the output. 
# 		the_good_stuff = the_good_stuff.replace(':', ',')
# 		final_line = str(str(chainID) + ',' + str(res_num) + ',' + str(the_good_stuff) + '\n')
# 		ramalyze_results.write(final_line)

# 	print((2 *'\n'), 'Ramachandran angle results written to', out_csv_name ,(2 * '\n'))
>>>>>>> ad95adbefe0997511532561a48002d6b2517df26

	# Generate chainID columns
	chain_IDs = [polypep.id] * len(chain_phis)

	# Add all polypep information to model DataFrame
	chain_summary = {
	'chainID' : chain_IDs,
	'residueName' : chain_resnames,
	'residueIndex' : chain_resindices,
	'phi' : chain_phis,
	'psi' : chain_psis,
	'type': chain_types
	}

	chain_summary = pd.DataFrame.from_dict(chain_summary)

	return chain_summary



def ModelDihedrals(model, model_num, iter_chains=True, chain_id=None):
	'''
	==============================================================================================
	Generates a Pandas DataFrame of phi/psi angles (and other information) from a given PDB model 
	==============================================================================================
	'''

	model_summaryDF = pd.DataFrame(columns=['chainID','residueName','residueIndex','phi','psi','type'])	

	# Iterate over all chains in model
	if iter_chains:

		for chain in model:
			chain_summaryDF = ChainSummary(chain)
			model_summaryDF = model_summaryDF.append(chain_summaryDF, ignore_index=True)

	# If multiple chains are present in the structure, default is to calculate dihedrals from all chains
	else:
		chain = model[chain_id]
		chain_summaryDF = ChainSummary(chain)
		model_summaryDF = model_summaryDF.append(chain_summaryDF, ignore_index=True)

<<<<<<< HEAD
	# Append model number information to final DataFrame
	model_ID_list = [model_num] * len(model_summaryDF)
	model_summaryDF.insert(loc=0, column='ModelID', value=model_ID_list)

	return model_summaryDF


=======
# CONDITIONAL CHECK TO ENSURE USER ENTERS THEIR PDB FILE NAME
# if len(sys.argv) == 1:
# 	print('\n Enter PDB file name: \n \n python3 CalculateDihedrals.py <file-name> \n')
# 	sys.exit()

# pdb_file_name = sys.argv[1]

# rama_angles_list = RamalyzeCommand(pdb_file_name)

# RamalyzeParseToCSV(rama_angles_list)
>>>>>>> ad95adbefe0997511532561a48002d6b2517df26

def ExtractDihedrals(pdb_file_name=None, iter_models=True, model_number=0, iter_chains=True, chain_id=None):
	'''
	==============================================================================================
	Generates a Pandas DataFrame of phi/psi angles (and other information) from a given PDB file 
	==============================================================================================
	'''
	inv_model = False

	# User parsed in a PDB file name
	if pdb_file_name != None:

		# Attempts to extract information from PDB file
		try:
			pdb_code = pdb_file_name[:-4]
			pdb_summaryDF = pd.DataFrame(columns=['ModelID','chainID','residueName','residueIndex','phi','psi','type'])

			# User did not parse in specific model: Iterate over all models in PDB object
			if iter_models:

				models = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_file_name)

				for model in models:

					model_dihedrals = ModelDihedrals(model, model_number, iter_chains, chain_id)
					pdb_summaryDF = pdb_summaryDF.append(model_dihedrals, ignore_index=True)

					model_number += 1

			# Specific model number parsed in by user
			else:
				try:
					model = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_file_name)[model_number]
					model_dihedrals = ModelDihedrals(model, model_number)
					pdb_summaryDF = pdb_summaryDF.append(model_dihedrals, ignore_index=True)
				# Invalid model number given 
				except:
					inv_model = True
					print('\n  ERROR: Invalid model number entered \n')
					exit()

			# Append PDB code information to final DataFrame
			pdb_list = [pdb_code] * len(pdb_summaryDF)
			pdb_summaryDF.insert(loc=0, column='PDBCode', value=pdb_list)

			return pdb_summaryDF

		# Invalid PDB file name given
		except:
			if not inv_model:
				print('\n  ERROR: Invalid PDB file \n ' )
				exit()
			else:
				exit()

	# No file name given
	else:
		print('\n ERROR: No PDB file specified \n')
		print(' Specify PDB file using: \n \n     --PDB /path_to_file/<filename.pdb> \n \n')
		exit()