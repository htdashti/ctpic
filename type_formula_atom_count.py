from aux_functions import carbohydrate_min_length, get_mass_of_a_path, get_atom_dic_of_apath
from aux_functions import add_carbohydrate_output


def check_formula_atom_count(formula_dic, atoms):
    output = []
    """
    monosaccharides:
     Parent monosaccharides are
        polyhydroxy aldehydes   H-[CHOH]n-CHO
        polyhydroxy ketones     H-[CHOH]n-CO-[CHOH]m-H
        with three or more carbon atoms.
    """
    if 'C' not in formula_dic or formula_dic['C'] < carbohydrate_min_length:
        return output
    if 'C' in formula_dic and 'O' in formula_dic and 'H' in formula_dic:
        # exact match
        if len(formula_dic) == 3:
            if formula_dic['C'] == formula_dic['O'] and formula_dic['H'] == 2*formula_dic['C']:
                atom_list = [_ for _ in range(len(atoms))]
                mass = get_mass_of_a_path(atom_list, atoms)
                full_path_atom_dic = get_atom_dic_of_apath(atom_list, atoms)
                output = add_carbohydrate_output(output, 'Atom count: polyhydroxy [CHOH]n', mass, formula_dic['C'], formula_dic['O'], atom_list, full_path_atom_dic)
    return output
