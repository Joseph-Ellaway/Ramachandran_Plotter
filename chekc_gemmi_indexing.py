


import gemmi

if __name__ == "__main__":

    ramas = {
		aa: [] for aa in [
			'LEU', 'ALA', 'GLY', 'VAL', 'GLU', 'SER', 'LYS', 'ASP', 'THR', 'ILE',
			'ARG', 'PRO', 'ASN', 'PHE', 'GLN', 'TYR', 'HIS', 'MET', 'CYS', 'TRP'
		]
	}


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

    # for chain in model.subchains():
    #     print(chain.subchain_id())
    #     print(dir(chain))

    # Show all chains in the model

    # chain_from_model = model['P']
#     chain_from_model = model.get_subchain('A')
# #
#     print(dir(chain_from_model))
#     print(chain_from_model.auth_seq_id_to_label)

    # for res in chain_from_model:
    #     pass

        # print(dir(res))
        # print(res.sifts_unp)
        # print(res.subchain)
        # print(res.seqid)
        # prev_res = chain_from_model.previous_residue(res)
        # next_res = chain_from_model.next_residue(res)

        # if prev_res and next_res and next_res.name != 'PRO':
        #     v = gemmi.calculate_phi_psi(prev_res, res, next_res)
        #     try:
        #         ramas[res.name].append(v)
        #     except KeyError:
                # pass
