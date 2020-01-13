from aux_functions import is_in_a_cycle, carbohydrate_chain_length_including_start_C, add_carbohydrate_output, \
    get_dic_of_number_of_atoms_in_path, get_mass_of_a_path, parse_atom_nghs, get_atom_dic_of_apath
import networkx as nx
import sys


def report(rep):
    if False:
        print(rep)


def find_c4(atoms):
    """
    CR4
    :param atoms:
    :return:
    """
    c4_group = []
    stop_ngh_trace = {}
    for c1_index in range(len(atoms)):
        if atoms[c1_index].get_atom_name() != "C":
            continue
        nghs_c1 = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1, ["H", "O"])
        if nums_c1["O"] == 0 and nums_c1["H"] == 0:
            c4_group.append(c1_index)
            stop_ngh_trace[c1_index] = []
            for _ in nghs_c1:
                stop_ngh_trace[c1_index].append(_[0])
    return c4_group, stop_ngh_trace


def find_cho_ccc(atoms):
    """
             C
            /
    R-C----C-H
      ||    \
      O     C
    or
      H      C
      |     /
    R-C----C-H
      |    \
      O     C
    or
      H      S
      |     /
    R-C----C-H
      |    \
      O     S
      example:
      (5R)-5-C-Cyclohexyl-5-C-phenyl-D-xylose
      2-Carb-16.4. Terminal substitution
      https://www.qmul.ac.uk/sbcs/iupac/2carb/16.html
    """
    cho_ccc_group = []
    stop_ngh_trace = {}
    for c1_index in range(len(atoms)):
        if atoms[c1_index].get_atom_name() != "C":
            continue
        nghs_c1 = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1, ["H", "O", "C"])
        if nums_c1["O"] == 1:
            for c2_index in nghs_list_c1["C"]:
                nghs_c2 = atoms[c2_index].get_ngh()
                nums_c2, nghs_list_c2 = parse_atom_nghs(nghs_c2, ["H", "O", "C"])
                if nums_c2["C"] >= 1 and nums_c2["H"] == 1 and nums_c2["O"] == 1:
                    cho_ccc_group.append(c1_index)
                    stop_ngh_trace[c1_index] = c2_index
    return cho_ccc_group, stop_ngh_trace


def find_r_ch_ch_r(atoms):
    """
    R-CH=CH-R
    :param atoms:
    :return:
    """
    ch_ch_group = []
    stop_ngh_trace = {}
    for c1_index in range(len(atoms)):
        if atoms[c1_index].get_atom_name() != "C":
            continue
        nghs_c1 = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1, ["H", "O", "C"])
        if len(nghs_c1) == 3:
            for c2_index in nghs_list_c1["C"]:
                if len(atoms[c2_index].get_ngh()) == 3:
                    ch_ch_group.append(c1_index)
                    stop_ngh_trace[c1_index] = c2_index
    return ch_ch_group, stop_ngh_trace


def find_cho(atoms):
    """
      H
      |
    R-C=O
    """
    cho_group = []
    for c1_index in range(len(atoms)):
        if atoms[c1_index].get_atom_name() != "C":
            continue
        nghs_c1 = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1, ["H", "O"])
        if len(nghs_c1) == 3 and nums_c1["H"] == 1 and nums_c1["O"] == 1:
            cho_group.append(c1_index)
    return cho_group


def find_chx_ox_r_groups(atoms):
    """
    H2
    |
    C-O-R
    |
    R
    or
    O
    |
    C-O-R
    |
    R
    :param atoms:
    :return:
    """
    chx_o_r_group = []
    for c1_index in range(len(atoms)):
        if atoms[c1_index].get_atom_name() != "C":
            continue
        nghs_c1 = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1, ["C", "H", "O"])
        if nums_c1["C"] == 1 and nums_c1["O"] != 0 and nums_c1["H"] >= 1:
            chx_o_r_group.append(c1_index)
        if nums_c1["C"] == 1 and nums_c1["O"] == 2:
            chx_o_r_group.append(c1_index)
    return chx_o_r_group


def find_ch3_groups(atoms):
    """
    C-H3
    |
    R
    :param atoms:
    :return:
    """
    ch3_group = []
    for c1_index in range(len(atoms)):
        if atoms[c1_index].get_atom_name() != "C":
            continue
        nghs_c1 = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1, ["C", "H", "O"])
        if nums_c1["H"] == 3:
            ch3_group.append(c1_index)
    return ch3_group


def find_conhoh_groups(atoms):
    """
    R--C---N----O
      ||   |    |
      O    H    H
    :param atoms:
    :return:
    """
    conhoh_groups = []
    for c1_index in range(len(atoms)):
        if atoms[c1_index].get_atom_name() != "C":
            continue
        nghs_c1 = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1, ["C", "H", "O", "N"])
        if len(nghs_c1) == 3 and nums_c1["O"] == 1 and nums_c1["N"] == 1:
            n1_index = nghs_list_c1["N"][0]
            nghs_n1 = atoms[n1_index].get_ngh()
            nums_n1, nghs_list_n1 = parse_atom_nghs(nghs_n1, ["C", "H", "O", "N"])
            if nums_n1["O"] == 1:
                conhoh_groups.append(c1_index)
    return conhoh_groups


def find_coch2oh_groups(atoms):
    """
    R--C---C----O
      ||   |    |
      O    H2   H
    :param atoms:
    :return:
    """
    conhoh_groups = []
    for c1_index in range(len(atoms)):
        if atoms[c1_index].get_atom_name() != "C":
            continue
        nghs_c1 = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1, ["C", "H", "O", "N"])
        if len(nghs_c1) == 3 and nums_c1["O"] == 1 and nums_c1["C"] == 2:
            for c2_index in nghs_list_c1["C"]:
                nghs_c2 = atoms[c2_index].get_ngh()
                nums_c2, nghs_list_c2 = parse_atom_nghs(nghs_c2, ["C", "H", "O", "N"])
                if nums_c2["O"] == 1 and nums_c2["H"] == 2:
                    o_index = nghs_list_c2["O"][0]
                    nghs_o = atoms[o_index].get_ngh()
                    nums_o, nghs_list_o = parse_atom_nghs(nghs_o, ["C", "H", "O", "N"])
                    if nums_o["H"] == 1 and nums_o["C"] == 1:
                        conhoh_groups.append(c1_index)
    return conhoh_groups


def find_ch2(atoms):
    """
    r-ch2
    example:1,2,4,5-Tetradeoxy-D-arabino-octa-1,4-dienitol
    :param atoms:
    :return:
    """
    ch2_group = []
    for c1_index in range(len(atoms)):
        if atoms[c1_index].get_atom_name() != "C":
            continue
        nghs_c1 = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1, ["C", "H", "O", "N"])
        if len(nghs_c1) == 3 and nums_c1["H"] == 2:
            ch2_group.append(c1_index)
    return ch2_group


def type_linear_structures_generic(molx, cycles, atoms):
    def add_subgroup(terminal_carbons, group):
        for _ in group:
            if _ not in terminal_carbons:
                terminal_carbons.append(_)
        return terminal_carbons

    def check_path_doesnt_have_ring(molx, cycles, start_index, end_index):
        paths = nx.all_simple_paths(molx, start_index, end_index)
        paths_dic = {}
        index = 0
        for a_path in paths:
            paths_dic[index] = {"len": len(a_path), "path": a_path}
            index += 1
        sorted(paths_dic.items(), key=lambda x: x[0])
        for a_path_dic in paths_dic:
            a_path = paths_dic[a_path_dic]["path"]
            if all([not is_in_a_cycle(cycles, _) for _ in a_path]):
                chain_atom_dic = get_dic_of_number_of_atoms_in_path(a_path, atoms)
                return [True, chain_atom_dic, a_path]
        return [False, {}, []]

    def count_num_o_in_a_path(a_path, atoms):
        count = 0
        for _ in a_path:
            nghs_c3 = atoms[_].get_ngh()
            nums_c3, nghs_list_c3 = parse_atom_nghs(nghs_c3, ['H', 'C', 'O', 'N'])
            if nums_c3['O'] != 0 or nums_c3['N'] != 0:
                count += 1
        return count
    output = []
    terminal_carbons = []
    terminal_carbons = add_subgroup(terminal_carbons, find_cho(atoms))
    terminal_carbons = add_subgroup(terminal_carbons, find_chx_ox_r_groups(atoms))
    terminal_carbons = add_subgroup(terminal_carbons, find_ch3_groups(atoms))
    terminal_carbons = add_subgroup(terminal_carbons, find_conhoh_groups(atoms))
    terminal_carbons = add_subgroup(terminal_carbons, find_coch2oh_groups(atoms))
    terminal_carbons = add_subgroup(terminal_carbons, find_ch2(atoms))
    cho_ccc_group, stop_ngh_trace_cho_ccc = find_cho_ccc(atoms)
    terminal_carbons = add_subgroup(terminal_carbons, cho_ccc_group)

    ch_ch_group, stop_ngh_trace_ch_ch = find_r_ch_ch_r(atoms)
    terminal_carbons = add_subgroup(terminal_carbons, ch_ch_group)

    c4_group, stop_ngh_trace_c4 = find_c4(atoms)
    terminal_carbons = add_subgroup(terminal_carbons, c4_group)

    terminal_carbons_outside_of_rings = []
    for _ in terminal_carbons:
        if not is_in_a_cycle(cycles, _):
            terminal_carbons_outside_of_rings.append(_)
    terminal_carbons = terminal_carbons_outside_of_rings

    report("terminal_carbons")
    report(terminal_carbons)
    for atom_index in terminal_carbons:
        # find paths to other terminals
        for another_terminal in terminal_carbons:
            if another_terminal == atom_index:
                continue
            [flag, chain_atom_dic, a_path] = check_path_doesnt_have_ring(molx, cycles, atom_index, another_terminal)
            num_o_connected_c = count_num_o_in_a_path(a_path, atoms)
            acceptable_number_of_missing_o = num_o_connected_c >= .8 * (len(a_path)-2)  # len path - head/tail terminal
            report("<<<<<<<<<<<<<<<<<<<<")
            report(flag)
            report(chain_atom_dic)
            report(a_path)
            report([num_o_connected_c, .8 * (len(a_path)-2)])  # len path - head/tail terminal

            if flag and chain_atom_dic['C'] >= carbohydrate_chain_length_including_start_C and acceptable_number_of_missing_o:
                path_mass = get_mass_of_a_path(a_path, atoms)
                nghs_to_be_ignored = []
                if atom_index in cho_ccc_group:
                    nghs_to_be_ignored.append(stop_ngh_trace_cho_ccc[atom_index])
                if another_terminal in cho_ccc_group:
                    nghs_to_be_ignored.append(stop_ngh_trace_cho_ccc[another_terminal])
                if atom_index in ch_ch_group:
                    nghs_to_be_ignored.append(stop_ngh_trace_ch_ch[atom_index])
                if another_terminal in ch_ch_group:
                    nghs_to_be_ignored.append(stop_ngh_trace_ch_ch[another_terminal])
                if atom_index in c4_group:
                    for _ in stop_ngh_trace_c4[atom_index]:
                        nghs_to_be_ignored.append(_)
                if another_terminal in c4_group:
                    for _ in stop_ngh_trace_c4[another_terminal]:
                        nghs_to_be_ignored.append(_)

                report("nghs_to_be_ignored:")
                report(nghs_to_be_ignored)
                full_path_atom_dic = get_atom_dic_of_apath(a_path, atoms, nghs_to_be_ignored)
                report(full_path_atom_dic)
                output = add_carbohydrate_output(output, 'generic saccharide chain',
                                                 path_mass, chain_atom_dic['C'], chain_atom_dic['O'], a_path,
                                                 full_path_atom_dic)

    return output