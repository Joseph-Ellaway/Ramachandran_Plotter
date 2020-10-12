
''' ================================================================================================================

		Deploy in the same directory as the .pdb from which the phenix.ramalyze results wish to be generated. 

	================================================================================================================ '''


	
import subprocess
from subprocess import PIPE
import sys



########################################################
#						FUNCTIONS 					   #
########################################################
def RamalyzeCommand(pdb_file):
	'''
	=======================================================
	Runs phenix.ramalyze command and processes the results as a list of strings. 
	=======================================================
	'''
	pdb_ramalyze_stdout = subprocess.run(['phenix.ramalyze', pdb_file], stdout=PIPE)
	pdb_ramalyze_stdout_list = pdb_ramalyze_stdout.stdout.decode().split('\n')
	rama_angles_list = pdb_ramalyze_stdout_list[1:-4]		# Removes the SUMMARY lines from the list
	return rama_angles_list




def RamalyzeParseToCSV(stdout_list):
	'''
	=======================================================
	Takes every element from the penix.ramalyse results list and parse them into a .csv file to be 
	read with Pandas in future graph making scripts.
	=======================================================
	'''
	out_csv_name = str(sys.argv[1][:-4] + '_DihedralAngles.csv')
	ramalyze_results = open(out_csv_name, "w")
	ramalyze_results.write(str('chainID,residue,residueID,score%,phi,psi,evaluation,type' + '\n'))		# Adds the PDB file name to the start of the description row. 
	ramalyze_results_listolists = []

	for residue_info in stdout_list:
		residue_info_split = list(residue_info.split(' '))
		ramalyze_results_listolists.append(residue_info_split)

	ramalyze_results_listolists_no_spaces = []

	for embedded_list in ramalyze_results_listolists:
		rama_info_no_spaces = []
		for item in embedded_list:
			if item != '':
				rama_info_no_spaces.append(item)
			else:
				pass

		ramalyze_results_listolists_no_spaces.append(rama_info_no_spaces)

	for rama_info_list_no_spaces in ramalyze_results_listolists_no_spaces:
		chainID = rama_info_list_no_spaces[0]
		res_num = rama_info_list_no_spaces[1]
		the_good_stuff = max(rama_info_list_no_spaces, key=len)		# len key used so to avoid getting high values, such as 'v's, in the output. 
		the_good_stuff = the_good_stuff.replace(':', ',')
		final_line = str(str(chainID) + ',' + str(res_num) + ',' + str(the_good_stuff) + '\n')
		ramalyze_results.write(final_line)

	print((2 *'\n'), 'Ramachandran angle results written to', out_csv_name ,(2 * '\n'))





########################################################
#				PROCESSING USER'S DATA				   #
########################################################

# CONDITIONAL CHECK TO ENSURE USER ENTERS THEIR PDB FILE NAME
if len(sys.argv) == 1:
	print('\n Enter PDB file name: \n \n python3 CalculateDihedrals.py <file-name> \n')
	sys.exit()

pdb_file_name = sys.argv[1]

rama_angles_list = RamalyzeCommand(pdb_file_name)

RamalyzeParseToCSV(rama_angles_list)

