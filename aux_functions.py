carbohydrate_min_length = 3
carbohydrate_chain_length_including_start_C = 3


def is_glycan(a_category):
    for _ in a_category["types"]:
        if "Glycan" in _:
            return True
    return False


def add_carbohydrate_output(output, type_, mass, num_c_atoms, num_o_atoms, atom_list, full_atom_dic):
    output.append({'type': type_,
                   'mass': mass,
                   'CO_count': [num_c_atoms, num_o_atoms],
                   'atom_indices': atom_list,
                   'full_atom_dic': full_atom_dic
                   })
    return output


def are_in_one_cycle(cycles, indices):  # indices: an array of indices
    for a_cycle in cycles:
        all_in = True
        for an_index in indices:
            if an_index not in a_cycle:
                all_in = False
        if all_in:
            return a_cycle
    return []


def is_in_a_cycle(cycles, index):
    for a_cycle in cycles:
        if index in a_cycle:
            return True
    return False


def there_is_double_bond(mol, atom_index1, atom_index2):
    atom_index1 += 1
    atom_index2 += 1
    bonds = mol.get_bonds()
    for a_bond in bonds:
        from_index, to_index, num_bonds, val1 = a_bond.get_info()
        if (atom_index1 == from_index and atom_index2 == to_index) or (atom_index1 == to_index and atom_index2 == from_index):
            if num_bonds == 2:
                return True
            else:
                return False


def parse_atom_nghs(nghs, atoms_of_interest):
    nums = {}
    nghs_list = {}
    for an_atom in atoms_of_interest:
        nums[an_atom] = 0
        nghs_list[an_atom] = []
    for a_ngh in nghs:
        if a_ngh[1] in nums:
            nums[a_ngh[1]] += 1
            nghs_list[a_ngh[1]].append(a_ngh[0])
    return nums, nghs_list


def get_dic_of_number_of_atoms_in_path(a_path, atoms):
    dic = {}
    for an_atom_index in a_path:
        curr_atom_name = atoms[an_atom_index].get_atom_name()
        if curr_atom_name in dic:
            dic[curr_atom_name] += 1
        else:
            dic[curr_atom_name] = 1
        curr_nghs = atoms[an_atom_index].get_ngh()
        cleaned_curr_nghs = []
        for a_pair in curr_nghs:
            if a_pair[0] not in a_path:
                cleaned_curr_nghs.append(a_pair)
        nghs_names_to_consider = ['H', 'C', 'O', 'N']
        nums, nghs_list = parse_atom_nghs(cleaned_curr_nghs, nghs_names_to_consider)
        """ we want to include num N's when counting O's """
        nums['O'] += nums['N']  # number of atoms
        for an_atom_name in nghs_names_to_consider:
            if an_atom_name == 'N':
                continue
            if an_atom_name not in dic:
                dic[an_atom_name] = nums[an_atom_name]
            else:
                dic[an_atom_name] += nums[an_atom_name]
    return dic


def get_mass_of_a_path(a_path, atoms):
    total_mass = 0
    for _ in range(len(a_path)):
        index = a_path[_]
        total_mass += atoms[index].get_mass()
        nghs_temp = atoms[index].get_ngh()
        nghs = [_ for _ in nghs_temp if _[0] not in a_path]
        nums, nghs_list = parse_atom_nghs(nghs, ['H', 'C', 'O', 'N'])
        # adding protons attached to a heavy atom in the path
        total_mass += sum([atoms[_].get_mass() for _ in nghs_list['H']])
        # adding mass of nghs and their protons
        for a_heavy_atom_type in ['C', 'O', 'N']:
            for _ in nghs_list[a_heavy_atom_type]:
                # mass of heavy atom
                total_mass += atoms[_].get_mass()
                # mass of its protons that are not in the path
                nghs_temp = atoms[_].get_ngh()
                nums, nghs_list = parse_atom_nghs(nghs_temp, ['H', 'C', 'O', 'N'])
                total_mass += sum([atoms[_].get_mass() for _ in nghs_list['H']])
    return total_mass


def get_atom_dic_of_apath(a_path, atoms, nghs_to_be_ignored):
    """ in this function, we get types and number of atoms in a path.
    we also continue every branch of these atoms until we reach to a non_ O, C, H atom
    """

    def find_acceptable_nghs(atoms, atom_index, a_path):
        """ returns nghs C or O atoms that are not in the a_path """
        if atoms[atom_index].get_atom_name() not in ['C', 'O']:  # we do not branch out from N
            return []
        nghs_temp = atoms[atom_index].get_ngh()
        nghs = [_ for _ in nghs_temp if _[0] not in a_path and _[0] not in nghs_to_be_ignored]
        return nghs

    def process_nghs(dic, seen_atoms, atoms, atom_index, a_path, ngh_level):
        ngh_level += 1
        if ngh_level > 3:  # and atoms[atom_index].get_atom_name() != 'O':
            return [dic, seen_atoms]
        nghs = find_acceptable_nghs(atoms, atom_index, a_path)
        if not nghs:
            return [dic, seen_atoms]
        if ngh_level > 2:
            current_acceptable_atom_names = ['H']
        else:
            current_acceptable_atom_names = acceptable_atom_names
        nums, nghs_list = parse_atom_nghs(nghs, current_acceptable_atom_names)
        if ngh_level >= 1 and atoms[atom_index].get_atom_name() == 'O':
            curr_nghs_list = {'H': nghs_list['H']}
        else:
            curr_nghs_list = nghs_list
        for _ in current_acceptable_atom_names:
            if _ not in curr_nghs_list:
                continue
            for __ in curr_nghs_list[_]:
                if __ not in seen_atoms and __ not in a_path:
                    if _ not in dic:
                        dic[_] = 0
                    dic[_] += 1
                    seen_atoms.append(__)
            for a_ngh_index in curr_nghs_list[_]:
                [dic, seen_atoms] = process_nghs(dic, seen_atoms, atoms, a_ngh_index, a_path, ngh_level)
        return [dic, seen_atoms]

    acceptable_atom_names = ['H', 'C', 'O', 'N']
    dic = {}
    seen_atoms = []
    for atom_index in a_path:
        if atom_index in seen_atoms:
            continue
        seen_atoms.append(atom_index)
        atom_name = atoms[atom_index].get_atom_name()  # atom can be one of the acceptable atom names
        if atom_name not in acceptable_atom_names:
            continue
        if atom_name not in dic:
            dic[atom_name] = 1
        else:
            dic[atom_name] += 1
        ngh_level = 0

        [dic, seen_atoms] = process_nghs(dic, seen_atoms, atoms, atom_index, a_path, ngh_level)
    dic['atoms_indices'] = seen_atoms  # [_+1 for _ in seen_atoms]
    return dic
