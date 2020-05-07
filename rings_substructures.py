from aux_functions import are_in_one_cycle, get_dic_of_number_of_atoms_in_path, \
    parse_atom_nghs, there_is_double_bond
import sys


def get_cor_groups(atoms, cycles):
    cor_groups = []
    interesting_atom_types = ["S", "O", "N", "F"]

    for atom_index in range(len(atoms)):
        if atoms[atom_index].get_atom_name() == "C":
            nghs_o_list = atoms[atom_index].get_ngh()
            nums, nghs_list = parse_atom_nghs(nghs_o_list, ['H', 'C', 'O', "F", "S", "N"])
            for __ in interesting_atom_types:
                for _ in nghs_list[__]:
                    if not are_in_one_cycle(cycles, [atom_index, _]) and atom_index not in cor_groups:
                        cor_groups.append(atom_index)
    return cor_groups


def get_c_ch2_or(atoms, cycles):
    """
    R
    |
    c-ch2-O-r
    |
    R
    :param atoms:
    :param cycles:
    :return:
    """
    c_ch2_or_group = []
    for c1_index in range(len(atoms)):
        if atoms[c1_index].get_atom_name() != "C":
            continue
        nghs = atoms[c1_index].get_ngh()
        nums, nghs_list = parse_atom_nghs(nghs, ['H', 'C', 'O'])
        for c2_index in nghs_list["C"]:
            if not are_in_one_cycle(cycles, [c1_index, c2_index]):
                nghs_c2 = atoms[c2_index].get_ngh()
                nums_c2, nghs_list_c2 = parse_atom_nghs(nghs_c2, ['H', 'C', 'O'])
                if nums_c2["O"] >= 1 and nums_c2["H"] == 2:
                    c_ch2_or_group.append(c1_index)
    return c_ch2_or_group


def get_c_chx_or(atoms, a_circle):
    """
    R
    |
    c-chx-O-r
    |
    R
    :param atoms:
    :param cycles:
    :return:
    """
    c_ch2_or_group = []
    for c1_index in range(len(atoms)):
        if atoms[c1_index].get_atom_name() != "C":
            continue
        nghs = atoms[c1_index].get_ngh()
        nums, nghs_list = parse_atom_nghs(nghs, ['H', 'C', 'O'])
        for c2_index in nghs_list["C"]:
            if c2_index not in a_circle:
                nghs_c2 = atoms[c2_index].get_ngh()
                nums_c2, nghs_list_c2 = parse_atom_nghs(nghs_c2, ['H', 'C', 'O'])
                if nums_c2["O"] >= 1 and nums_c2["H"] == 1:
                    c_ch2_or_group.append(c1_index)
    return c_ch2_or_group


def find_co_double_bonds(atoms, mol):
    """ looking for
        R
        |
        C=O
        |
        R
    """
    CO_list = []
    for _ in range(len(atoms)):
        name = atoms[_].get_atom_name()
        if name != 'C':
            continue
        nghs_c = atoms[_].get_ngh()
        nums_c, nghs_list_c = parse_atom_nghs(nghs_c, ['H', 'O'])
        if nums_c["O"] == 1 and nums_c["H"] == 0:
            if there_is_double_bond(mol, nghs_list_c["O"][0], _):
                CO_list.append(_)
    return CO_list


def find_C_CH3(atoms):
    """
    want to find
    R
    |
    C-CH3
    |
    R
    """
    CCH3_list = []
    for _ in range(len(atoms)):
        name = atoms[_].get_atom_name()
        if name != 'C':
            continue
        nghs_c = atoms[_].get_ngh()
        nums_c, nghs_list_c = parse_atom_nghs(nghs_c, ["C", "N"])
        for __ in nghs_list_c["C"]:
            curr_nghs_c = atoms[__].get_ngh()
            curr_nums_c, curr_nghs_list_c = parse_atom_nghs(curr_nghs_c, ['H'])
            if curr_nums_c["H"] == 3:
                CCH3_list.append(_)

    return CCH3_list


def find_terminal_rings(mol, atoms, cycles):
    """ example: 4'a-carbathymidine
        2-Carb-34.2. Replacement by carbon
        https://www.qmul.ac.uk/sbcs/iupac/2carb/34.html
    """
    CCH3_list = find_C_CH3(atoms)
    CO_list = find_co_double_bonds(atoms, mol)
    def is_ring_thymidine(a_cycle):
        ring_dic = get_dic_of_number_of_atoms_in_path(a_cycle, atoms)
        if "N" in ring_dic and ring_dic["N"] == 2:
            num_co_in_ring = 0
            num_c_ch3_in_ring = 0
            for _ in a_cycle:
                name = atoms[_].get_atom_name()
                if name == 'C' and _ in CO_list:
                    num_co_in_ring += 1
                if name == "C" and _ in CCH3_list:
                    num_c_ch3_in_ring += 1
            if num_co_in_ring == 2 and num_c_ch3_in_ring:
                return True
        return False
    c_with_terminal_ring = []
    for _ in range(len(atoms)):
        name = atoms[_].get_atom_name()
        if name != 'C':
            continue
        nghs_c = atoms[_].get_ngh()
        nums_c, nghs_list_c = parse_atom_nghs(nghs_c, ["N"])
        if nums_c["N"] == 1:
            N_index = nghs_list_c["N"][0]
            for a_cycle in cycles:
                if N_index in a_cycle and is_ring_thymidine(a_cycle):
                    c_with_terminal_ring.append(_)
    return c_with_terminal_ring


def find_C_with_N_terminals(atoms):
    """
    want to find
    R
    |
    C-NH3
    |
    R
    or
    R
    |
    C-NH2
    |
    R
    or
    R
    |
    C-NH-CO-CH3
    |
    R
    """
    CNH3_list = []
    for _ in range(len(atoms)):
        name = atoms[_].get_atom_name()
        if name != 'C':
            continue
        nghs_c = atoms[_].get_ngh()
        nums_c, nghs_list_c = parse_atom_nghs(nghs_c, ["N"])
        for __ in nghs_list_c["N"]:
            curr_nghs_c = atoms[__].get_ngh()
            curr_nums_c, curr_nghs_list_c = parse_atom_nghs(curr_nghs_c, ['H', "C"])
            if curr_nums_c["C"] == 1 and (curr_nums_c["H"] == 3 or curr_nums_c["H"] == 2):
                CNH3_list.append(_)
            if curr_nums_c["C"] == 2:
                for another_c in curr_nghs_list_c["C"]:
                    if another_c == _:
                        continue
                    another_nghs_c = atoms[another_c].get_ngh()
                    another_nums_c, another_nghs_list_c = parse_atom_nghs(another_nghs_c, ['H', "C", "O"])
                    if another_nums_c["O"] == 1 and another_nums_c["H"] == 0 and another_nums_c["C"] == 1:
                        third_nghs_c = atoms[another_nghs_list_c["C"][0]].get_ngh()
                        third_nums_c, third_nghs_list_c = parse_atom_nghs(third_nghs_c, ['H', "C", "O"])
                        if third_nums_c["H"] == 3:
                            CNH3_list.append(_)

    return CNH3_list


def find_c_ch2_o_p_o(atoms, cycles):
    """
    looking for
    R       R
    |       |
    C-CH2-O-P=O
    |       |
    R       R
    :param atoms:
    :return:
    """
    c_ch2_o_p_o = []
    for c1_index in range(len(atoms)):
        if atoms[c1_index].get_atom_name() != "C":
            continue
        nghs_c1 = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1, ['C'])
        for c2_index in nghs_list_c1["C"]:
            if not are_in_one_cycle(cycles, [c1_index, c2_index]):
                nghs_c2 = atoms[c2_index].get_ngh()
                nums_c2, nghs_list_c2 = parse_atom_nghs(nghs_c2, ['C', "H", "O"])
                if nums_c2["C"] == 1 and nums_c2["O"] == 1:
                    o1_index = nghs_list_c2["O"][0]
                    nghs_o1 = atoms[o1_index].get_ngh()
                    nums_o1, nghs_list_o1 = parse_atom_nghs(nghs_o1, ["P"])
                    if nums_o1["P"] == 1:
                        p1_index = nghs_list_o1["P"][0]
                        nghs_p1 = atoms[p1_index].get_ngh()
                        nums_p1, nghs_list_p1 = parse_atom_nghs(nghs_p1, ["O"])
                        if nums_p1["O"] >= 2:
                            c_ch2_o_p_o.append(c1_index)

    return c_ch2_o_p_o


def find_c_o(atoms, a_circle):
    """
          R
         /
    R-C-O
      |
      R
    :param atoms:
    :return:
    """
    c_o_group = []
    for c1_index in range(len(atoms)):
        if atoms[c1_index].get_atom_name() != "C":
            continue
        nghs_c1 = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1, ["O"])
        for o1_index in nghs_list_c1["O"]:
            if o1_index not in a_circle:
                c_o_group.append(c1_index)
    return c_o_group


def find_c_coor(atoms, a_circle, cycles):
    c_coor_group = []
    for c1_index in a_circle:
        if atoms[c1_index].get_atom_name() != "C":
            continue
        nghs_c1 = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1, ["C"])
        for c2_index in nghs_list_c1["C"]:
            in_circle = False
            for a_cycle in cycles:
                if c2_index in a_cycle:
                    in_circle = True
            if not in_circle:
                nghs_c2 = atoms[c2_index].get_ngh()
                nums_c2, nghs_list_c2 = parse_atom_nghs(nghs_c2, ["O"])
                if nums_c2["O"] >= 1:
                    c_coor_group.append(c1_index)
    return c_coor_group


def find_xylos_rings(mol, atoms, cycles):  #, coh_groups, CH2OH_groups, CH2OR_groups, CH2OH_CO_groups, COOH_R_in_chain, COXx_R_groups, carbon_with_terminal_rings, CNH3, cor_groups, c_ch2_o_p_o_groups):
    # an example is X4S
    def report(rep):
        if False:
            print(rep)

    def add_terminal_groups(a_set, group):
        output = [_ for _ in a_set]
        for _ in group:
            if _ not in a_set:
                output.append(_)
        return output
    terminal_subgroups = []
    terminal_subgroups = add_terminal_groups(terminal_subgroups, get_c_ch2_or(atoms, cycles))
    terminal_subgroups = add_terminal_groups(terminal_subgroups, get_cor_groups(atoms, cycles))
    terminal_subgroups = add_terminal_groups(terminal_subgroups, find_terminal_rings(mol, atoms, cycles))
    terminal_subgroups = add_terminal_groups(terminal_subgroups, find_C_with_N_terminals(atoms))
    terminal_subgroups = add_terminal_groups(terminal_subgroups, find_c_ch2_o_p_o(atoms, cycles))
    report("terminal_subgroups")
    report(terminal_subgroups)
    xylos_groups = []
    messages_xylos_groups = []
    for a_circle in cycles:
        messages = []
        if len(a_circle) >= 7:
            continue
        curr_terminal_subgroups = add_terminal_groups(terminal_subgroups, get_c_chx_or(atoms, a_circle))
        curr_terminal_subgroups = add_terminal_groups(curr_terminal_subgroups, find_c_coor(atoms, a_circle, cycles))
        curr_terminal_subgroups = add_terminal_groups(curr_terminal_subgroups, find_c_o(atoms, a_circle))
        report("<<<<")
        report([_ for _ in a_circle])
        report([_ for _ in curr_terminal_subgroups])
        report(terminal_subgroups)
        report(">>>>")
        # does the circle have an O
        has_an_O = False
        ring_o_list = []
        ring_other_list = []
        for atom_name in [atoms[_].get_atom_name() for _ in a_circle]:
            if atom_name == "O":
                ring_o_list.append(atom_name)
            if atom_name in ["S", "P", "N", "Se"]:
                ring_other_list.append(atom_name)
        if ring_o_list:
            has_an_O = True
        else:
            if ring_other_list:
                c_count = 0
                num_four_bound_carbons = 0
                for _ in a_circle:
                    if atoms[_].get_atom_name() == "C":
                        c_count += 1
                        nghs_list = atoms[_].get_ngh()
                        if len(nghs_list) >= 4:
                            num_four_bound_carbons += 1
                if c_count == num_four_bound_carbons:
                    has_an_O = True
                    messages.append("ring O substitution to 'S', 'P', 'N', or 'Se'")
        # for ring O subs to C, all of the atoms in the circle should have 4 bonds. Otherwise cyclitols?!
        if not has_an_O:
            num_four_bound_carbons = 0
            for _ in a_circle:
                if atoms[_].get_atom_name() == "C":
                    nghs_list = atoms[_].get_ngh()
                    if len(nghs_list) >= 4:
                        num_four_bound_carbons += 1
            if num_four_bound_carbons == len(a_circle):
                has_an_O = True
                messages.append("ring O substitution to 'C'")
        report("has O: %d" % has_an_O)
        if has_an_O:
            num_coh = 0
            c_counter = 0
            for an_atom in a_circle:
                if atoms[an_atom].get_atom_name() == "C":
                    c_counter += 1
                if an_atom in curr_terminal_subgroups:
                    num_coh += 1
            report("num_coh: %d" % num_coh)
            if num_coh >= 3:
                if c_counter < 3:
                    print(a_circle)
                    print("c_counter < 3")
                    sys.exit(3)
                xylos_groups.append(a_circle)
    return xylos_groups