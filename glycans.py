from aux_functions import parse_atom_nghs
from aux_functions import add_carbohydrate_output
from aux_functions import get_dic_of_number_of_atoms_in_path
from aux_functions import get_mass_of_a_path
from aux_functions import get_atom_dic_of_apath
from aux_functions import is_in_a_cycle


def get_glycan_atom_info(glycan_atoms, atoms):
    cycle_atoms_dic = get_dic_of_number_of_atoms_in_path(glycan_atoms, atoms)
    cycle_mass = get_mass_of_a_path(glycan_atoms, atoms)
    full_path_atom_dic = get_atom_dic_of_apath(glycan_atoms, atoms, [])
    return cycle_atoms_dic, cycle_mass, full_path_atom_dic


def type_ring_glycan_fucose(cycles, atoms, R_CHOR_R_fragments, R_CHOR_R_fragments_atom_lists,
                            R_CHOH_R_fragments, R_CHOH_R_fragments_atom_lists,
                            Nindiex_R_NH_CO_R_in_chain, Nindiex_R_NH_CO_R_in_chain_atoms_list,
                            CH2OH_groups, CH2OH_groups_atoms_list):
    output = []
    if not cycles:
        return output
    for a_cycle in cycles:
        cycle_atoms = []
        for _ in a_cycle:
            cycle_atoms.append([_, atoms[_].get_atom_name()])
        nums, nghs = parse_atom_nghs(cycle_atoms, ['H', 'C', 'O'])
        if nums["O"] != 1 or nums["C"] != 5:
            continue
        in_r_choh_r_list = 0
        num_NH_CO_R = 0
        num_CH2OH = 0
        num_CH2_in_ring = 0
        num_CHOH_CHOH_side_chain = 0
        additional_atoms = {"R_CHOH_R_fragments": [],
                            "R_CHOR_R_fragments": [],
                            "Nindiex_R_NH_CO_R": [],
                            "CH2OH_groups": [],
                            "R_CHOH_CHOH_R": []}

        def add_to_additional_atoms(array, dic_term, array_flag):
            for _ in array:
                if array_flag:
                    __ = _[0]
                else:
                    __ = _
                if __ not in additional_atoms[dic_term]:
                    additional_atoms[dic_term].append(__)

        for c_index in nghs["C"]:
            if c_index in R_CHOH_R_fragments or c_index in R_CHOR_R_fragments:
                in_r_choh_r_list += 1
                if c_index in R_CHOH_R_fragments:
                    add_to_additional_atoms(R_CHOH_R_fragments_atom_lists[c_index], "R_CHOH_R_fragments", array_flag=False)
                if c_index in R_CHOR_R_fragments:
                    add_to_additional_atoms(R_CHOR_R_fragments_atom_lists[c_index], "R_CHOR_R_fragments", array_flag=False)
            else:
                tmp_c = atoms[c_index].get_ngh()
                tmp_nums, tmp_nghs = parse_atom_nghs(tmp_c, ['H', 'C', 'O', 'N'])
                if tmp_nums["H"] == 1:
                    """_R_NH_CO_R_"""
                    for _ in tmp_nghs["N"]:
                        if _ in Nindiex_R_NH_CO_R_in_chain:
                            num_NH_CO_R += 1
                            add_to_additional_atoms(Nindiex_R_NH_CO_R_in_chain_atoms_list[_], "Nindiex_R_NH_CO_R", array_flag=False)
                    """CH2OH"""
                    for _ in tmp_nghs["C"]:
                        if _ in CH2OH_groups:
                            num_CH2OH += 1
                            add_to_additional_atoms(CH2OH_groups_atoms_list[_], "CH2OH_groups", array_flag=False)

                    """R_CHOH_CHOH_R"""
                    for c2_index in tmp_nghs["C"]:
                        if not is_in_a_cycle([a_cycle], c2_index):
                            c2_nghs_ngh_list = atoms[c2_index].get_ngh()
                            c2_nums, c2_nghs = parse_atom_nghs(c2_nghs_ngh_list, ['H', 'C', 'O', 'N'])
                            if c2_nums["H"] == 1 and c2_nums["O"] == 1 and c2_nums["C"] > 1:
                                for c3_index in c2_nghs["C"]:
                                    if not is_in_a_cycle([a_cycle], c3_index):
                                        c3_nghs_ngh_list = atoms[c3_index].get_ngh()
                                        c3_nums, c3_nghs = parse_atom_nghs(c3_nghs_ngh_list, ['H', 'C', 'O', 'N'])
                                        if c3_nums["H"] >= 1 and c3_nums["O"] >= 1:
                                            num_CHOH_CHOH_side_chain += 1
                                            add_to_additional_atoms(c2_nghs_ngh_list, "R_CHOH_CHOH_R", array_flag=True)
                                            add_to_additional_atoms(c3_nghs_ngh_list, "R_CHOH_CHOH_R", array_flag=True)
                if tmp_nums["H"] == 2:
                    num_CH2_in_ring += 1
        if num_CH2_in_ring > 1:
            # has more than 1 CH2 in the ring. Cannot be a glycan
            continue

        glycan_atoms = a_cycle
        for tags in ["Nindiex_R_NH_CO_R", "R_CHOH_R_fragments"]:
            for _ in additional_atoms[tags]:
                if _ not in glycan_atoms:
                    glycan_atoms.append(_)

        if in_r_choh_r_list >= 4:  # "Is Hex; Fucose", gluctose, glucose, mannose, N-acetylgalactosamine
            cycle_atoms_dic, cycle_mass, full_path_atom_dic = get_glycan_atom_info(glycan_atoms, atoms)
            output = add_carbohydrate_output(output, "Glycan; Hex", cycle_mass,
                                             cycle_atoms_dic["C"], cycle_atoms_dic["O"],
                                             glycan_atoms, full_path_atom_dic)

        if in_r_choh_r_list == 3 and (num_NH_CO_R == 1 and num_CH2OH == 1):  # could be N-acetylehexosamine
            cycle_atoms_dic, cycle_mass, full_path_atom_dic = get_glycan_atom_info(glycan_atoms, atoms)
            output = add_carbohydrate_output(output, "Glycan; N-acetylehexosamine", cycle_mass,
                                             cycle_atoms_dic["C"], cycle_atoms_dic["O"],
                                             glycan_atoms, full_path_atom_dic)
        if in_r_choh_r_list == 1 and \
            num_CH2_in_ring == 1 and \
                num_NH_CO_R == 1 and \
                num_CHOH_CHOH_side_chain >= 1:  # could be N-acetyleneuraminic_acid
            cycle_atoms_dic, cycle_mass, full_path_atom_dic = get_glycan_atom_info(glycan_atoms, atoms)
            output = add_carbohydrate_output(output, "Glycan; N-acetyleneuraminic_acid", cycle_mass,
                                             cycle_atoms_dic["C"], cycle_atoms_dic["O"],
                                             glycan_atoms, full_path_atom_dic)
    return output


def type_ring_glycan_xylose(xylos_rings, coh_groups, atoms, mol):
    output = []
    for a_xyl_ring in xylos_rings:
        cycle_atoms_dic, cycle_mass, full_path_atom_dic = get_glycan_atom_info(a_xyl_ring, atoms)
        output = add_carbohydrate_output(output, "Glycan; xylo derivatives", cycle_mass,
                                         cycle_atoms_dic["C"], cycle_atoms_dic["O"],
                                         a_xyl_ring, full_path_atom_dic)
    return output