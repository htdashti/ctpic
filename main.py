import os
import sys
import compound
import networkx as nx
import json
import calculate_penalties
import requests

from unsaturated_monosaccharides import unsaturated_carb_17
from linear_substructures import type_linear_structures_generic
from rings_substructures import find_xylos_rings

from find_fragments import find_CH2OH_in_chain
from find_fragments import find_R_CHOH_R_in_ring
from find_fragments import find_R_NH_CO_R_in_chain
from find_fragments import find_R_CHOR_R_in_ring
from find_fragments import get_coh_groups

from glycans import type_ring_glycan_fucose
from glycans import type_ring_glycan_xylose


def group_carbohydrates_by_atom_indices(input_set_of_carbohydrates):
    groups_of_carbohydrates = []
    seen = [False] * len(input_set_of_carbohydrates)
    for _ in range(len(input_set_of_carbohydrates)):
        if not seen[_]:
            seen[_] = True
            new_group = {'type': [input_set_of_carbohydrates[_]['type']],
                         'mass': input_set_of_carbohydrates[_]['mass'],
                         'CO_count': input_set_of_carbohydrates[_]['CO_count'],
                         'atom_indices': input_set_of_carbohydrates[_]['atom_indices'],
                         'full_atom_dic': input_set_of_carbohydrates[_]['full_atom_dic']
                         }
            for __ in range(len(input_set_of_carbohydrates)):
                if not seen[__]:
                    if set(input_set_of_carbohydrates[_]['atom_indices']) == \
                            set(input_set_of_carbohydrates[__]['atom_indices']) or \
                            set(input_set_of_carbohydrates[_]['full_atom_dic']['atoms_indices']) == \
                            set(input_set_of_carbohydrates[__]['full_atom_dic']['atoms_indices']):
                        if input_set_of_carbohydrates[__]['type'] not in new_group['type']:
                            new_group['type'].append(input_set_of_carbohydrates[__]['type'])
                        seen[__] = True
            groups_of_carbohydrates.append(new_group)
    return groups_of_carbohydrates


def convert_table_to_dic(mol, formula_str, formula_dic, classes, sdf_file_path):
    dic = {'mol_mass': mol.get_mass(), 'mol_num_atoms': len(mol.get_atoms()),
           'Formula': formula_str, 'formula_dic': formula_dic,
           'Carbohydrate_categories': [], 'fname': sdf_file_path}
    grouped_input_set = []
    for a_row in classes:
        for a_category in a_row[1]:
            grouped_input_set.append(a_category)
    groups_of_carbohydrates = group_carbohydrates_by_atom_indices(grouped_input_set)
    for a_category in groups_of_carbohydrates:
        cat_dic = {'types': []}
        for a_type in a_category['type']:
            cat_dic['types'].append(a_type)
        cat_dic['mass'] = a_category['mass']
        cat_dic['full_atom_dic'] = a_category['full_atom_dic']
        cat_dic['num_atoms_in_carbohydrate'] = len(a_category['full_atom_dic']['atoms_indices'])
        cat_dic['atom_indices'] = [_ + 1 for _ in a_category['atom_indices']]
        dic['Carbohydrate_categories'].append(cat_dic)
    return dic


def extract_carbohydrate_flagments(mol, sdf_file_path):
    atoms = mol.get_atoms()
    molx = mol.get_networkx()
    cycles = list(nx.cycle_basis(molx))

    formula_dic = mol.get_formula_dic()
    formula_str = mol.get_formula_str().replace(' ', '')

    coh_groups = get_coh_groups(mol, atoms)
    CH2OH_groups, CH2OH_groups_atoms_list = find_CH2OH_in_chain(atoms, cycles)

    xylos_rings = find_xylos_rings(mol, atoms, cycles)

    Nindiex_R_NH_CO_R_in_chain, Nindiex_R_NH_CO_R_in_chain_atoms_list = find_R_NH_CO_R_in_chain(mol, atoms, cycles)

    R_CHOH_R_fragments, R_CHOH_R_fragments_atom_lists = find_R_CHOH_R_in_ring(atoms, cycles)
    R_CHOR_R_fragments, R_CHOR_R_fragments_atom_lists = find_R_CHOR_R_in_ring(atoms, cycles)

    classes = [['type_ring_glycan_fucose', type_ring_glycan_fucose(cycles, atoms, R_CHOR_R_fragments,
                                                                   R_CHOR_R_fragments_atom_lists,
                                                                   R_CHOH_R_fragments,
                                                                   R_CHOH_R_fragments_atom_lists,
                                                                   Nindiex_R_NH_CO_R_in_chain,
                                                                   Nindiex_R_NH_CO_R_in_chain_atoms_list,
                                                                   CH2OH_groups, CH2OH_groups_atoms_list)],
               ['type_ring_glycan_xylose', type_ring_glycan_xylose(xylos_rings, coh_groups, atoms, mol)],
               ['unsaturated', unsaturated_carb_17(mol, atoms, cycles)],
               ["linear generic", type_linear_structures_generic(molx, cycles, atoms)]]

    entries_carbohydrate_dic = convert_table_to_dic(mol, formula_str, formula_dic, classes, sdf_file_path)
    return entries_carbohydrate_dic


def calculate_compound_penalty(mol, init_formula_dic):
    curr_formula_dic = mol.get_formula_dic()
    num_heavy_atoms_cmp = 0
    for an_atom_type in curr_formula_dic:
        if an_atom_type != 'H':
            num_heavy_atoms_cmp += curr_formula_dic[an_atom_type]
    curr_num_heavy_atoms = num_heavy_atoms_cmp

    num_heavy_atoms_carb = 0
    for an_atom_type in init_formula_dic:
        if an_atom_type != 'H':
            num_heavy_atoms_carb += init_formula_dic[an_atom_type]
    init_num_heavy_atoms = num_heavy_atoms_carb
    cmp_penalty = curr_num_heavy_atoms / init_num_heavy_atoms
    return cmp_penalty


def run_on_a_structure(output_root, sdf_file_path):
    def report(rep):
        if False:
            print(rep)
    mol = compound.load_sdf(sdf_file_path, 'SDF')
    init_formula_dic = mol.get_formula_dic()

    entries_carbohydrate_dic = extract_carbohydrate_flagments(mol, sdf_file_path)
    report('Carbohydrate_categories')
    for _ in entries_carbohydrate_dic['Carbohydrate_categories']:
        report(_)

    if not entries_carbohydrate_dic['Carbohydrate_categories']:
        # done with possible carbohydrates
        if "C" not in init_formula_dic or "O" not in init_formula_dic:
            cmp_penalty = 1
        else:
            cmp_penalty = calculate_compound_penalty(mol, init_formula_dic)
        entries_carbohydrate_dic["best_penalty"] = 1
        entries_carbohydrate_dic["Carbohydrate_categories"] = []
        report("cmp_penalty: %f" % cmp_penalty)
    else:
        best_penalty, best_fragment, updated_new_categories = \
            calculate_penalties.process(entries_carbohydrate_dic, mol.get_atoms())
        to_be_removed_atom_indices = []
        report("updated_new_categories: ")
        for _ in updated_new_categories:
            report(_)
        for an_updated_cat in updated_new_categories:
            for __ in an_updated_cat["full_atom_dic"]["atoms_indices"]:
                to_be_removed_atom_indices.append(__)
        report("num atoms in carb fragments: %d" % len(to_be_removed_atom_indices))
        mol.remove_atoms(to_be_removed_atom_indices)
        mol.calculate_formula()
        cmp_penalty = calculate_compound_penalty(mol, init_formula_dic)
        report("cmp_penalty: %f" % cmp_penalty)
        report("best penalty: %f" % best_penalty)
        entries_carbohydrate_dic["best_penalty"] = best_penalty
        entries_carbohydrate_dic["Carbohydrate_categories"] = updated_new_categories

    entries_carbohydrate_dic["cmp_penalty"] = cmp_penalty
    json.dump(entries_carbohydrate_dic, open(os.path.join(output_root, "output.json"), "w"), indent=4)

    pdb_ligands_similar_cmp = get_most_similar_structures(sdf_file_path)

    os.system('rm -f *.pyc')
    os.chdir(output_root)
    json.dump(pdb_ligands_similar_cmp, open("similar_cmp_pdb_map.json", "w"), indent=4)
    os.system('zip -q output_web.zip similar_cmp_pdb_map.json output.json')
    return entries_carbohydrate_dic["best_penalty"]


def get_alatis_structure(sdf_file_path):
    req = requests.post('http://alatis.nmrfam.wisc.edu/upload',
                        data={'format': 'sdf', 'response_type': 'json', 'project_2_to_3': 'on', 'add_hydrogens': 'on'},
                        files={'infile': open(sdf_file_path, 'r')}).json()
    if req["status"] == "error":
        req["structure"] = open(sdf_file_path, 'r').read()
        req["html_url"] = "ALATIS could not process the input file"
        req["inchi"] = "ALATIS could not process the input file"
    return req["structure"], req["html_url"], req["inchi"], req["status"]


def get_most_similar_structures(alatis_output_file_name):
    ligans_to_pdb_links = json.load(open("ligands_links_to_pdb.json", "r"))

    def get_pdb_protein_ids(ligand_id):
        pdb_ids = []
        if ligand_id in ligans_to_pdb_links:
            pdb_ids = ligans_to_pdb_links[ligand_id]
        return pdb_ids
    command = "obabel pdb_ligands.fs -O obabel_similar_cmps.sdf -s%s -at5" % alatis_output_file_name
    os.system(command)
    content = open("obabel_similar_cmps.sdf", "r").read().split("$$$$\n")
    similar_entries = []
    for a_file in content:
        lines = a_file.split("\n")
        similar_entries.append(lines[0])
    pdb_map = {}
    for _ in similar_entries:
        if _:
            pdb_map[_] = get_pdb_protein_ids(_)
    return pdb_map


if __name__ == "__main__":
    sdf_file_path = sys.argv[1]
    alatis_structure, alatis_url, alatis_inchi, alatis_status = get_alatis_structure(sdf_file_path)
    alatis_output_file_name = sdf_file_path  # "alatis_called_to_generate_structure.sdf"
    fout = open(alatis_output_file_name, "w")
    fout.write(alatis_structure)
    fout.close()

    msg = {'sdf_file_path': alatis_output_file_name,
           "alatis_url": alatis_url,
           "alatis_inchi": alatis_inchi,
           "input_structure": sdf_file_path,
           "alatis_status": alatis_status,
           "error": "",
           "warning": ""}
    if alatis_status == "error":
        msg["warning"] = "ALATIS could not process the input file"
    try:
        best_penalty = run_on_a_structure("./", sdf_file_path)
        msg["run_message"] = "run_msg"
        try:
            similar_cmp_pdb_map = get_most_similar_structures(alatis_output_file_name)
            json.dump(similar_cmp_pdb_map, open("similar_cmp_pdb_map.json", "w"), indent=4)
        except:
            msg["error"] = "Could not find similar compounds"
    except:
        msg["error"] = "Something went wrong! Could not finish processing the input file"
    json.dump(msg, open("msg.json", "w"), indent=4)
    os.system('zip -q output_web.zip msg.json output.json similar_cmp_pdb_map.json %s' % sdf_file_path)