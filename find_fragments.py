from aux_functions import is_in_a_cycle, are_in_one_cycle, there_is_double_bond, parse_atom_nghs, get_dic_of_number_of_atoms_in_path


def get_coh_groups(mol, atoms):
    coh_groups = []
    for atom_index in range(len(atoms)):
        if atoms[atom_index].get_atom_name() == "O":
            nghs_o_list = atoms[atom_index].get_ngh()
            nums, nghs_list = parse_atom_nghs(nghs_o_list, ['H', 'C', 'O'])
            if nums["H"] == 1 and nums["C"] == 1:
                coh_groups.append(nghs_list["C"][0])
    return coh_groups


def get_cor_groups(atoms, cycles):
    cor_groups = []
    for atom_index in range(len(atoms)):
        if atoms[atom_index].get_atom_name() == "C":
            nghs_o_list = atoms[atom_index].get_ngh()
            nums, nghs_list = parse_atom_nghs(nghs_o_list, ['H', 'C', 'O', "F"])
            for _ in nghs_list["O"]:
                if not are_in_one_cycle(cycles, [atom_index, _]):
                    cor_groups.append(atom_index)
            for _ in nghs_list["F"]:
                if not are_in_one_cycle(cycles, [atom_index, _]) and atom_index not in cor_groups:
                    cor_groups.append(atom_index)
    return cor_groups



def find_R_CHOR_R_in_line(atoms, cycles):
    """
        looking for
        case 1-A:

              R
             /
            O
            '
         R--C(1)--R
            '
            H

        in line
    """
    # case 1-A
    CHOR_groups = []
    for c1_index in range(len(atoms)):
        name = atoms[c1_index].get_atom_name()
        if name != 'C':
            continue
        if not is_in_a_cycle(cycles, c1_index):
            continue
        nghs_c1_list = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1_list, ['H', 'C', 'O'])
        if nums_c1["H"] != 1:
            continue
        for o_index in nghs_list_c1["O"]:
            if not are_in_one_cycle(cycles, [c1_index, o_index]):  # not in a same cycle
                CHOR_groups.append(c1_index)
                break
    return CHOR_groups


def find_R_CH3_in_line(atoms, cycles):
    """
        looking for
        case 1-A:
            H
            '
         R--C(1)--H
            '
            H

        in line
    """
    # case 1-A
    CH3_groups = []
    for c1_index in range(len(atoms)):
        name = atoms[c1_index].get_atom_name()
        if name != 'C':
            continue
        if is_in_a_cycle(cycles, c1_index):
            continue
        nghs_c1_list = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1_list, ['H', 'C', 'O'])
        if nums_c1["H"] == 3:
            CH3_groups.append(c1_index)
    return CH3_groups


def find_R_CHOR_R_in_ring(atoms, cycles):
    """
        looking for
        case 1-A:

              R
             /
            O
            '
         R--C(1)--R
            '
            H

        in rings
    """
    # case 1-A
    ring_fragments_atom_lists = {}
    ring_fragments = []
    for c1_index in range(len(atoms)):
        name = atoms[c1_index].get_atom_name()
        if name != 'C':
            continue
        if not is_in_a_cycle(cycles, c1_index):
            continue
        nghs_c1_list = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1_list, ['H', 'C', 'O'])
        if nums_c1["H"] != 1:
            continue
        for o_index in nghs_list_c1["O"]:
            if not are_in_one_cycle(cycles, [c1_index, o_index]):  # not in a same cycle
                nghs_o_list = atoms[o_index].get_ngh()
                """ it should not be a double bond to C1 (i.e. len(nghs_o_list) > 1). However, because for composite 
                carbohydrates we iteratively remove the best fragment from a molecule, len(nghs_o_list) could be 1 """
                if True:  #len(nghs_o_list) > 1:  #
                    ring_fragments.append(c1_index)
                    ring_fragments_atom_lists[c1_index] = []
                    for _ in nghs_c1_list:
                        ring_fragments_atom_lists[c1_index].append(_[0])
                    for _ in nghs_o_list:
                        ring_fragments_atom_lists[c1_index].append(_[0])
    return ring_fragments, ring_fragments_atom_lists


def find_R_CHOH_R_in_ring(atoms, cycles):
    """
        looking for
        case 1-A:

            O-H
            '
         R--C(1)--R
            '
            H

        in rings
    """
    # case 1-A
    ring_fragments_atom_lists = {}
    ring_fragments = []
    for c1_index in range(len(atoms)):
        name = atoms[c1_index].get_atom_name()
        if name != 'C':
            continue
        if not is_in_a_cycle(cycles, c1_index):
            continue
        nghs_c1_list = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1_list, ['H', 'C', 'O'])
        if nums_c1["H"] != 1:
            continue
        for o_index in nghs_list_c1["O"]:
            if not are_in_one_cycle(cycles, [c1_index, o_index]):  # not in a same cycle
                nghs_o_list = atoms[o_index].get_ngh()
                nums_o, nghs_list_o = parse_atom_nghs(nghs_o_list, ['H', 'C', 'O'])
                if nums_o["C"] == 1 and nums_o["H"] == 1:
                    ring_fragments.append(c1_index)
                    ring_fragments_atom_lists[c1_index] = []
                    for _ in nghs_c1_list:
                        ring_fragments_atom_lists[c1_index].append(_[0])
                    for _ in nghs_o_list:
                        ring_fragments_atom_lists[c1_index].append(_[0])
    return ring_fragments, ring_fragments_atom_lists


def find_R_CH2_CR2_R_in_ring(atoms, cycles):
    """
        looking for
        case 1-A:
                   R1
            H      '
            '      '
        R---C(2)---C(1)---R2
            '      '
            H      '
                   R3
        in rings
    """
    # case 1-A
    ring_fragments = []
    for c2_index in range(len(atoms)):
        name = atoms[c2_index].get_atom_name()
        if name != 'C':
            continue
        nghs_c2 = atoms[c2_index].get_ngh()
        nums_c2, nghs_list_c2 = parse_atom_nghs(nghs_c2, ['H', 'C', 'O'])
        if nums_c2['H'] != 2:
            continue
        for c1_index in nghs_list_c2['C']:
            cycle_0 = are_in_one_cycle(cycles, [c2_index, c1_index])
            if cycle_0:
                ring_fragments.append([c2_index, c1_index, cycle_0])

    return ring_fragments


def find_R_N_CR2_R_in_ring(atoms, cycles):
    """
        looking for
        case 1-A:
            H      R1
            '      '
            '      '
        R---N(2)---C(1)---R2
                   '
                   '
                   R3
        or
        case 1-B:
            H      R1
            '      '
            '      '
        R---P(2)---C(1)---R2
                   '
                   '
                   R3
        in rings
    """
    # case 1-A
    ring_fragments = []
    for n2_index in range(len(atoms)):
        name = atoms[n2_index].get_atom_name()
        if name != 'N':
            continue
        nghs_n2 = atoms[n2_index].get_ngh()
        nums_n2, nghs_list_n2 = parse_atom_nghs(nghs_n2, ['H', 'C', 'O'])
        for c1_index in nghs_list_n2['C']:
            cycle_0 = are_in_one_cycle(cycles, [n2_index, c1_index])
            if cycle_0:
                ring_fragments.append([n2_index, c1_index, cycle_0])
    # case 1-B
    for p2_index in range(len(atoms)):
        name = atoms[p2_index].get_atom_name()
        if name != 'P':
            continue
        nghs_p2 = atoms[p2_index].get_ngh()
        nums_p2, nghs_list_p2 = parse_atom_nghs(nghs_p2, ['H', 'C', 'O'])
        if nums_p2['H'] != 1:
            continue
        for c1_index in nghs_list_p2['C']:
            cycle_0 = are_in_one_cycle(cycles, [p2_index, c1_index])
            if cycle_0:
                ring_fragments.append([p2_index, c1_index, cycle_0])
    return ring_fragments


def find_R_O_CR2_R_in_ring(atoms, cycles):
    """
        looking for
                   R1
                   '
                   '
        R---O(2)---C(1)---R2
                   '
                   '
                   R3
        in rings
    """
    ring_fragments = []
    for o2_index in range(len(atoms)):
        name = atoms[o2_index].get_atom_name()
        if name != 'O':
            continue
        nghs_o2 = atoms[o2_index].get_ngh()
        nums_o2, nghs_list_o2 = parse_atom_nghs(nghs_o2, ['H', 'C', 'O'])
        for c1_index in nghs_list_o2['C']:
            cycle_0 = are_in_one_cycle(cycles, [o2_index, c1_index])
            if cycle_0:
                ring_fragments.append([o2_index, c1_index, cycle_0])
    return ring_fragments


def find_CH2OR_in_chain(atoms, cycles):
    """ this function finds terminal CH2OR that C is not in a cycle
                   R2
                   '
                   O(6)
                   '   H
                   '  /
               R---C(5)---H
    """

    end_carbon_indices = []
    for _ in range(len(atoms)):
        name = atoms[_].get_atom_name()
        if name != 'C' or is_in_a_cycle(cycles, _):
            continue
        nghs_c5 = atoms[_].get_ngh()
        nums_c5, nghs_list_c5 = parse_atom_nghs(nghs_c5, ['H', 'C', 'O'])
        if nums_c5['H'] == 2 and nums_c5['O'] == 1:
            o6_index = nghs_list_c5['O'][0]
            nghs_o6 = atoms[o6_index].get_ngh()
            nums_o6, nghs_list_o6 = parse_atom_nghs(nghs_o6, ['H', 'C', 'O'])
            if len(nghs_o6) == 2 and nums_o6['H'] == 0:
                end_carbon_indices.append(_)
    return end_carbon_indices


def find_CH2OH_in_chain(atoms, cycles):
    """ this function finds terminal CH2OH that C is not in a cycle
                   H
                   '
                   O(6)
                   '   H
                   '  /
               R---C(5)---H
    """

    end_carbon_indices = []
    end_carbon_indices_atom_list = {}
    for _ in range(len(atoms)):
        name = atoms[_].get_atom_name()
        if name != 'C' or is_in_a_cycle(cycles, _):
            continue
        nghs_c5 = atoms[_].get_ngh()
        nums_c5, nghs_list_c5 = parse_atom_nghs(nghs_c5, ['H', 'C', 'O'])
        if nums_c5['H'] == 2 and nums_c5['O'] == 1:
            o6_index = nghs_list_c5['O'][0]
            nghs_o6 = atoms[o6_index].get_ngh()
            nums_o6, nghs_list_o6 = parse_atom_nghs(nghs_o6, ['H', 'C', 'O'])
            if len(nghs_o6) == 2 and nums_o6['H'] == 1 and nums_o6['C'] == 1:
                end_carbon_indices.append(_)
                end_carbon_indices_atom_list[_] = []
                for __ in nghs_c5:
                    ___ = __[0]
                    if ___ not in end_carbon_indices_atom_list[_]:
                        end_carbon_indices_atom_list[_].append(___)
                for __ in nghs_o6:
                    ___ = __[0]
                    if ___ not in end_carbon_indices_atom_list[_]:
                        end_carbon_indices_atom_list[_].append(___)
    return end_carbon_indices, end_carbon_indices_atom_list


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


def find_R_CH2_R_groups(atoms):
    """ this function finds terminal CH2OH that C is not in a cycle
                   H
                   '
               R---C(1)---R
                   '
                   H
    """

    end_carbon_indices = []
    for _ in range(len(atoms)):
        name = atoms[_].get_atom_name()
        if name != 'C':
            continue
        nghs_c1 = atoms[_].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1, ['H', 'C', 'O'])
        if nums_c1['H'] == 2:
                end_carbon_indices.append(_)
    return end_carbon_indices


def find_R_O_R_groups(atoms, cycles):
    """
               R---O(1)---R

               R is not H
               O1 is in a cycle
    """
    end_carbon_indices = []
    for atom_index in range(len(atoms)):
        name = atoms[atom_index].get_atom_name()
        if (name != 'O' and name != "S") or not is_in_a_cycle(cycles, atom_index):
            continue
        nghs_o1 = atoms[atom_index].get_ngh()
        nums_o1, nghs_list_o1 = parse_atom_nghs(nghs_o1, ['H', 'C', 'O'])
        if nums_o1['H'] == 0:
                end_carbon_indices.append(atom_index)
    return end_carbon_indices


def find_CHO_in_chain(atoms, cycles):
    """ this function finds terminal CH2OH that C is not in a cycle
                   O(4)
                   "
                   "
               R---C(3)---H
    """
    end_carbon_indices = []
    for _ in range(len(atoms)):
        name = atoms[_].get_atom_name()
        if name != 'C' or is_in_a_cycle(cycles, _):
            continue
        nghs_c3 = atoms[_].get_ngh()
        nums_c3, nghs_list_c3 = parse_atom_nghs(nghs_c3, ['H', 'C', 'O'])
        if nums_c3['H'] == 1 and nums_c3['O'] == 1:
            o4_index = nghs_list_c3['O'][0]
            nghs_o4 = atoms[o4_index].get_ngh()
            if len(nghs_o4) == 1:
                end_carbon_indices.append(_)
    return end_carbon_indices


def find_CH2OH_CO_R_in_chain(atoms, cycles):
    """ these structures shouldn't be in rings
                          H
                          '
                          O(2)   O(4)
                      H   '      "
                        \ '      "
                      H---C(1)---C(3)---
    """
    end_carbon_indices = []
    for c1_index in range(len(atoms)):
        name = atoms[c1_index].get_atom_name()
        if name != 'C' or is_in_a_cycle(cycles, c1_index):
            continue
        nghs_c1 = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1, ['H', 'C', 'O'])
        if len(nghs_c1) == 4 and nums_c1['H'] == 2 and nums_c1['O'] == 1 and nums_c1['C'] == 1:
            o2_is_ok = False
            o2_index = nghs_list_c1['O'][0]
            nghs_o2 = atoms[o2_index].get_ngh()
            nums_o2, nghs_list_o2 = parse_atom_nghs(nghs_o2, ['H', 'C', 'O'])
            if len(nghs_o2) == 2 and nums_o2['C'] == 1 and nums_o2['H'] == 1:
                o2_is_ok = True
            if not o2_is_ok:
                continue
            c3_index = nghs_list_c1['C'][0]
            nghs_c3 = atoms[c3_index].get_ngh()
            nums_c3, nghs_list_c3 = parse_atom_nghs(nghs_c3, ['H', 'C', 'O'])
            if len(nghs_c3) == 3 and nums_c3['O'] == 1 and nums_c3['H'] == 0:
                o4_index = nghs_list_c3['O'][0]
                nghs_o4 = atoms[o4_index].get_ngh()
                if len(nghs_o4) == 1:
                    end_carbon_indices.append([c1_index, c3_index])
    return end_carbon_indices


def find_R_NH_CO_R_in_chain(mol, atoms, cycles):
    """ these structures shouldn't be in rings
                            O(2)
                     H      "
                     '      "
                R--- N(0)---C(1)---R
    """
    nitrogen_indices_atoms_list = {}
    nitrogen_indices = []
    for N0_index in range(len(atoms)):
        name = atoms[N0_index].get_atom_name()
        if name != 'N' or is_in_a_cycle(cycles, N0_index):
            continue
        nghs_n1 = atoms[N0_index].get_ngh()
        nums_n1, nghs_list_n1 = parse_atom_nghs(nghs_n1, ['H', 'C', 'O'])
        if nums_n1["H"] != 1 or nums_n1["C"] == 0 or nums_n1["O"] != 0:  # confirmed R-NH-C1 part
            continue
        for c1_index in nghs_list_n1["C"]:
            nghs_c1 = atoms[c1_index].get_ngh()
            nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1, ['H', 'C', 'O', 'N'])
            num_nghs = 0
            for _ in nums_c1:
                num_nghs += nums_c1[_]
            if num_nghs == 3 and nums_c1['O'] == 1:
                if there_is_double_bond(mol, c1_index, nghs_list_c1['O'][0]):  # confirmed C=O and -CO-
                    nitrogen_indices.append(N0_index)
                    nitrogen_indices_atoms_list[N0_index] = []
                    for _ in nghs_n1:
                        __ = _[0]
                        if __ not in nitrogen_indices_atoms_list[N0_index]:
                            nitrogen_indices_atoms_list[N0_index].append(__)
                    for _ in nghs_c1:
                        __ = _[0]
                        if __ not in nitrogen_indices_atoms_list[N0_index]:
                            nitrogen_indices_atoms_list[N0_index].append(__)
    return nitrogen_indices, nitrogen_indices_atoms_list


def find_R_CO_R_in_chain(mol, atoms, cycles):
    """ these structures shouldn't be in rings
                         O(2)
                         "
                         "
                     R---C(1)---R
    """
    end_carbon_indices = []
    for c1_index in range(len(atoms)):
        name = atoms[c1_index].get_atom_name()
        if name != 'C' or is_in_a_cycle(cycles, c1_index):
            continue
        nghs_c1 = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1, ['H', 'C', 'O'])
        if len(nghs_c1) == 3 and nums_c1['O'] == 1:
            if there_is_double_bond(mol, c1_index, nghs_list_c1['O'][0]):
                end_carbon_indices.append(c1_index)

    return end_carbon_indices


def find_chain_C_where_CO_in_ring(atoms, cycles):
    """
              O(2)--
              '     '
              '     '
          R---C(1)--R
              '
              '
              H
    """
    end_carbon_indices = []
    for c1_index in range(len(atoms)):
        name = atoms[c1_index].get_atom_name()
        if name != 'C' or not is_in_a_cycle(cycles, c1_index):
            continue
        nghs_c1 = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1, ['H', 'C', 'O'])
        if nums_c1['O'] == 1 and nums_c1['H'] == 1:
            if are_in_one_cycle(cycles, [c1_index, nghs_list_c1['O'][0]]):
                end_carbon_indices.append(c1_index)
    return end_carbon_indices


def find_COOH_R_in_chain(mol, atoms, cycles):
    """
              H
              '
              O(2)
              '
              '
       O(0)===C(1)--R
    """
    end_carbon_indices = []
    for c1_index in range(len(atoms)):
        name = atoms[c1_index].get_atom_name()
        if name != 'C' or is_in_a_cycle(cycles, c1_index):
            continue
        nghs_c1 = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1, ['H', 'C', 'O'])
        if len(nghs_c1) == 3 and nums_c1['O'] == 2:
            o0_index = -1
            o2_index = -1
            for an_o_index in nghs_list_c1['O']:
                if there_is_double_bond(mol, an_o_index, c1_index):
                    o0_index = an_o_index
                else:
                    nghs_o2 = atoms[an_o_index].get_ngh()
                    nums_o2, nghs_list_o2 = parse_atom_nghs(nghs_o2, ['H', 'C', 'O'])
                    if len(nghs_o2) == 2 and nums_o2['C'] == 1 and nums_o2['H'] == 1:
                        o2_index = an_o_index
            if o0_index != -1 and o2_index != -1:
                end_carbon_indices.append(c1_index)
    return end_carbon_indices


def find_COXx_R_in_chain(mol, atoms, cycles):
    """
    case 1:
              x
              '
              O(2)
              '
              '
       O(0)===C(1)--R

    case 2:
              x
              '
              N(2)
              '
              '
       O(0)===C(1)--R
    """
    end_carbon_indices = []
    for c1_index in range(len(atoms)):
        name = atoms[c1_index].get_atom_name()
        if name != 'C' or is_in_a_cycle(cycles, c1_index):
            continue
        nghs_c1 = atoms[c1_index].get_ngh()
        nums_c1, nghs_list_c1 = parse_atom_nghs(nghs_c1, ['H', 'C', 'O', 'N'])

        # case 1
        if len(nghs_c1) == 3 and nums_c1['O'] == 2:
            o0_index = -1
            o2_index = -1
            for an_o_index in nghs_list_c1['O']:
                if there_is_double_bond(mol, an_o_index, c1_index):
                    o0_index = an_o_index
                else:
                    nghs_o2 = atoms[an_o_index].get_ngh()
                    nums_o2, nghs_list_o2 = parse_atom_nghs(nghs_o2, ['H', 'C', 'O'])
                    if len(nghs_o2) == 2 and nums_o2['C'] == 1 and nums_o2['H'] == 0:  # x should not be H
                        o2_index = an_o_index
            if o0_index != -1 and o2_index != -1:
                end_carbon_indices.append(c1_index)
        # case 2
        if len(nghs_c1) == 3 and nums_c1['O'] == 1 and nums_c1['N'] == 1:
            o0_is_ok = False
            o0_index = nghs_list_c1['O'][0]
            if there_is_double_bond(mol, o0_index, c1_index):
                o0_is_ok = True

            n2_is_ok = False
            n2_index = nghs_list_c1['N'][0]
            nghs_n2 = atoms[n2_index].get_ngh()
            nums_n2, nghs_list_n2 = parse_atom_nghs(nghs_n2, ['H', 'C', 'O'])
            if nums_n2['C'] == 1:
                n2_is_ok = True
            if o0_is_ok and n2_is_ok:
                end_carbon_indices.append(c1_index)
    return end_carbon_indices