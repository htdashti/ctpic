from aux_functions import is_glycan


def get_num_atoms(full_atom_dic, atom_name):
    num = 0
    if atom_name in full_atom_dic:
        num = full_atom_dic[atom_name]
    return num


def calculate_carb_penalty_of_a_category(carb_atom_dic):
    heavy_atom_penalty = .1
    less_important_heavy_atom_penalty = .02  # this is for when num O is higher than number of C
    proton_penanlty = .05
    less_important_proton_penanlty = .01  # when |O| == |C|
    # carbohydrate penalty
    num_carbons = carb_atom_dic['C']  # number of carbons
    num_oxygens = get_num_atoms(carb_atom_dic, 'O')
    num_hydrogen = get_num_atoms(carb_atom_dic, 'H')
    if num_oxygens == num_carbons:
        carb_penalty = abs(2 * num_carbons - num_hydrogen) * less_important_proton_penanlty
    else:
        if num_oxygens > num_carbons:  # this case is better than |O| < |C|
            carb_penalty = abs(num_carbons - num_oxygens) * less_important_heavy_atom_penalty
            carb_penalty += min(
                [abs(2 * num_carbons - num_hydrogen), abs(2 * num_oxygens - num_hydrogen)]) * proton_penanlty
        else:
            carb_penalty = abs(num_carbons - num_oxygens) * heavy_atom_penalty
            carb_penalty += abs(2 * num_carbons - num_hydrogen) * proton_penanlty
    max_carb_penalty = max([num_carbons, num_oxygens]) * heavy_atom_penalty  # heavy atoms
    max_carb_penalty += max([2 * num_carbons, 2 * num_oxygens]) * proton_penanlty
    """ for non carbohydrate atoms we add a pentalty, for every appearance of the atoms """

    additional_atom_types = [_ for _ in carb_atom_dic if _ != 'C' and _ != 'O' and _ != 'H' and _ != 'atoms_indices']
    for an_atom in additional_atom_types:
        carb_penalty += carb_atom_dic[an_atom] * heavy_atom_penalty
        max_carb_penalty += carb_atom_dic[an_atom] * heavy_atom_penalty
    if max_carb_penalty == 0:  # there is no C, N, O, H in the fragment
        carb_penalty = 1
    else:
        carb_penalty = min([carb_penalty / max_carb_penalty, 1])
    return carb_penalty


def update_carb_atom_dic(carb_atom_dic, best_fragment_atom_dic, atoms):
    new_atom_indices = []
    for _ in carb_atom_dic["atoms_indices"]:
        if _ not in best_fragment_atom_dic["atoms_indices"]:
            new_atom_indices.append(_)

    new_dic = {"C": 0, "N": 0, "H": 0, "O": 0}
    for _ in new_atom_indices:
        atom_name = atoms[_].get_atom_name()
        if atom_name not in new_dic:
            new_dic[atom_name] = 0
        new_dic[atom_name] += 1
    new_dic["atoms_indices"] = new_atom_indices
    return new_dic


def process(data, atoms):
    best_penalty = 1
    best_fragment = {}
    to_be_removed_glycan_atoms = {'atoms_indices': []}

    def add_to_glycan_atoms(carb_atom_dic):
        # print(carb_atom_dic)
        for _ in carb_atom_dic['atoms_indices']:
            if _ not in to_be_removed_glycan_atoms['atoms_indices']:
                to_be_removed_glycan_atoms['atoms_indices'].append(_)
        for _ in to_be_removed_glycan_atoms['atoms_indices']:
            atom_name = atoms[_].get_atom_name()
            if atom_name not in to_be_removed_glycan_atoms:
                to_be_removed_glycan_atoms[atom_name] = 1
            else:
                to_be_removed_glycan_atoms[atom_name] += 1

    def is_the_best_category(a_category, full_atom_dic):
        for _ in a_category["full_atom_dic"]["atoms_indices"]:
            if _ not in full_atom_dic["atoms_indices"]:
                return False
        for _ in full_atom_dic["atoms_indices"]:
            if _ not in a_category["full_atom_dic"]["atoms_indices"]:
                return False
        return True
    for __ in range(len(data['Carbohydrate_categories'])):
        # for a_category in data['Carbohydrate_categories']:
        a_category = data['Carbohydrate_categories'][__]
        # this contain the number of atoms in a path that is identified as a carbohydrate
        carb_atom_dic = a_category['full_atom_dic']
        if is_glycan(a_category):
            add_to_glycan_atoms(carb_atom_dic)
            carb_penalty = 0
        else:
            carb_penalty = calculate_carb_penalty_of_a_category(carb_atom_dic)
        data['Carbohydrate_categories'][__]["fragment_penalty"] = carb_penalty
        if carb_penalty <= best_penalty:
            best_penalty = carb_penalty
            best_fragment = a_category
    new_categories = []

    for a_category in data['Carbohydrate_categories']:
        if is_glycan(a_category):
            a_category["cmp_penalty"] = 0
            new_categories.append(a_category)
            continue
        if is_the_best_category(a_category, best_fragment['full_atom_dic']):
            a_category["cmp_penalty"] = best_penalty
            new_categories.append(a_category)
            continue
        carb_atom_dic = a_category['full_atom_dic']
        updated_carb_atom_dic = update_carb_atom_dic(carb_atom_dic, best_fragment['full_atom_dic'], atoms)
        updated_carb_atom_dic = update_carb_atom_dic(updated_carb_atom_dic, to_be_removed_glycan_atoms, atoms)
        if not updated_carb_atom_dic["atoms_indices"]:  # all of its atoms have been removed by the best or glycans
            continue
        updated_penalty = calculate_carb_penalty_of_a_category(updated_carb_atom_dic)
        a_category["atoms_indices"] = updated_carb_atom_dic
        a_category["cmp_penalty"] = updated_penalty
        if updated_penalty < .7:
            new_categories.append(a_category)
    return best_penalty, best_fragment, new_categories
