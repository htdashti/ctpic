from aux_functions import is_in_a_cycle, add_carbohydrate_output, get_mass_of_a_path, parse_atom_nghs, get_atom_dic_of_apath


def unsaturated_carb_17(mol, atoms, cycles):
    output = []
    """
    case: 2-Carb-17. Unsaturated monosaccharides
        C(1)-C(2)-O(3)
        one of the C's should have a double or triple bonds
    """
    unsaturated_atoms_list = []
    bonds = mol.get_bonds()
    for a_bond in bonds:
        [from_index, to_index, num_bonds, val1] = a_bond.get_info()
        from_index -= 1
        to_index -= 1
        if num_bonds > 1 and atoms[from_index].get_atom_name() != 'O' and atoms[to_index].get_atom_name() != 'O':
            unsaturated_atoms_list.append(from_index)
            unsaturated_atoms_list.append(to_index)
    if not unsaturated_atoms_list:
        # there is no unsaturated atom
        return output
    for c1_index in unsaturated_atoms_list:
        name = atoms[c1_index].get_atom_name()
        if name != 'C' or is_in_a_cycle(cycles, c1_index):
            continue
        nghs_c1 = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1, ['H', 'C', 'O'])
        for c2_index in nghs_list_c1['C']:
            nghs_c2 = atoms[c2_index].get_ngh()
            nums_c2, nghs_list_c2 = parse_atom_nghs(nghs_c2, ['H', 'C', 'O'])
            for o3_index in nghs_list_c2['O']:
                nghs_o3 = atoms[o3_index].get_ngh()
                nums_o3, nghs_list_o3 = parse_atom_nghs(nghs_c2, ['H', 'C', 'O'])
                if len(nghs_o3) == 2 and nums_o3['C'] == 2:
                    atom_list = [c1_index, c2_index, o3_index]
                    mass = get_mass_of_a_path(atom_list, atoms)
                    full_path_atom_dic = get_atom_dic_of_apath(atom_list, atoms, [])
                    output = add_carbohydrate_output(output, '2-Carb-17. Unsaturated monosaccharides', mass, 2, 1, atom_list, full_path_atom_dic)
    return output
