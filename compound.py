import os
import known_elements
import networkx as nx

report = True


class atom:
    def __init__(self, name):
        self.name = name
        self.x = 0
        self.y = 0
        self.z = 0
        self.charge = 0
        self.mass = 0
        self.num_attached_proton = 0
        self.nghs = []

    def set_coord(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def add_ngh(self, ngh_index, ngh_name):
        self.nghs.append([ngh_index, ngh_name])

    def get_ngh(self):
        return self.nghs

    def set_charge(self, charge):
        self.charge = charge

    def set_mass(self, mass):
        self.mass = mass

    def get_atom_name(self):
        return self.name

    def get_mass(self):
        return self.mass

    def get_info(self):
        return [self.name, self.x, self.y, self.z, self.charge, self.mass]

    def get_num_attached_protons(self):
        return self.num_attached_proton

    def set_num_attached_protons(self, val):
        self.num_attached_proton = val

    def get_mass_properties(self):
        return self.mass, self.num_attached_proton, self.name


class bond:
    def __init__(self, from_index, to_index, num_bonds):
        self.from_index = from_index  # int
        self.to_index = to_index  # int
        self.num_bonds = num_bonds
        self.val1 = 0  # reserved for kcal/mol/r^2

    def set_val1(self, val1):
        self.val1 = val1

    def get_from_to_indices(self):
        return self.from_index, self.to_index

    def get_info(self):
        return [self.from_index, self.to_index, self.num_bonds, self.val1]


class molecule:
    def __init__(self):
        self.atoms = []
        self.bonds = []
        self.formula_dic = {}
        self.formula_str = ''
        self.mass = 0
        self.graph = []
        self.num_heavy_atoms = 0
        self.min_non_proton_mass = 100
        self.num_nodes_in_network = 0

    def remove_atoms(self, atom_indices_from_1):
        def get_reduced(x):
            smaller_indices_index = 0
            for _ in range(x-1):
                if _+1 in atom_indices_from_1:
                    smaller_indices_index += 1
            return smaller_indices_index

        current_bonds = self.bonds
        self.bonds = []
        for a_bond in current_bonds:
            x, y, num_bonds, val1 = a_bond.get_info()
            if not (x in atom_indices_from_1 or y in atom_indices_from_1):
                self.bonds.append(bond(x-get_reduced(x), y-get_reduced(y), num_bonds))
        current_atoms = self.atoms
        self.atoms = []
        for index in range(len(current_atoms)):
            if index+1 not in atom_indices_from_1:
                [name, x, y, z, _, _] = current_atoms[index].get_info()
                self.add_an_atom(name, x, y, z)

    def add_an_atom(self, atom_name, x, y, z):
        if atom_name != 'H':
            self.num_heavy_atoms += 1
        atom_num = len(self.atoms) + 1
        new_atom = atom(atom_name)
        new_atom.set_coord(x, y, z)
        self.atoms.append(new_atom)
        return atom_num

    def add_a_bond(self, from_index, to_index, num_bonds):
        bond_num = len(self.bonds) + 1
        self.bonds.append(bond(from_index, to_index, num_bonds))
        return bond_num

    def get_num_nodes_in_network(self):
        return self.num_nodes_in_network

    def get_min_non_proton_mass(self):
        return self.min_non_proton_mass

    def get_atoms(self):
        return self.atoms

    def get_bonds(self):
        return self.bonds

    def get_formula_str(self):
        return self.formula_str

    def get_formula_dic(self):
        return self.formula_dic

    def get_mass(self):
        return self.mass

    def get_networkx(self):
        return self.graph

    def calculate_num_protons_attached_to_heavy_atoms(self):
        for a_bond in self.bonds:
            [from_index, to_index, num_bonds, val1] = a_bond.get_info()
            from_index = from_index - 1
            to_index = to_index - 1
            if self.atoms[from_index].get_atom_name() != 'H' and self.atoms[to_index].get_atom_name() == 'H':
                self.atoms[from_index].set_num_attached_protons(self.atoms[from_index].get_num_attached_protons() + 1)
            if self.atoms[from_index].get_atom_name() == 'H' and self.atoms[to_index].get_atom_name() != 'H':
                self.atoms[to_index].set_num_attached_protons(self.atoms[to_index].get_num_attached_protons() + 1)

    def calculate_formula(self):
        formula_dic = {}
        atoms = self.get_atoms()
        for an_atom in atoms:
            name = an_atom.get_atom_name()
            if name in formula_dic:
                formula_dic[name] += 1
            else:
                formula_dic[name] = 1
        formula_str = ''
        for a in sorted(formula_dic):
            formula_str += a + str(formula_dic[a]) + ' '
        self.formula_dic = formula_dic
        self.formula_str = formula_str

    def calculate_networkx(self):
        graph = nx.Graph()
        nodes = []
        for i in range(len(self.atoms)):
            iter_atom_name = self.atoms[i].get_atom_name()
            if iter_atom_name != 'H':
                nodes.append(i)
        num_nodes_in_network = len(nodes)
        graph.add_nodes_from(nodes)
        edge_list = []
        for a_bond in self.bonds:
            [from_index, to_index, num_bonds, val1] = a_bond.get_info()
            from_index = from_index-1
            to_index = to_index-1
            if self.atoms[from_index].get_atom_name() != 'H' and self.atoms[to_index].get_atom_name() != 'H':
                edge_list.append((from_index, to_index))
        graph.add_edges_from(edge_list)
        self.graph = graph
        self.num_nodes_in_network = num_nodes_in_network

    def set_mass(self):
        total_mass = 0
        min_non_proton_mass = 1000
        el_mass = known_elements.get_elements_mass()
        for index in range(len(self.atoms)):
            a_name = self.atoms[index].get_atom_name()
            if a_name in el_mass:
                c_el_mass = el_mass[a_name]
                self.atoms[index].set_mass(c_el_mass)
                total_mass += c_el_mass
                if a_name != 'H' and c_el_mass < min_non_proton_mass:
                    min_non_proton_mass = c_el_mass
            else:
                print('****Error: atom mass for ' + a_name + ' ' + str(index) + ' was not found')
        self.mass = total_mass
        self.min_non_proton_mass = min_non_proton_mass

    def set_nghs(self):
        for a_bond in self.bonds:
            from_index, to_index, num_bonds, val1 = a_bond.get_info()
            from_index -= 1
            to_index -= 1
            self.atoms[from_index].add_ngh(to_index, self.atoms[to_index].get_atom_name())
            self.atoms[to_index].add_ngh(from_index, self.atoms[from_index].get_atom_name())

    def print_graph(self):
        print(self.graph.nodes())
        print(self.graph.edges())

    def print_general_info(self):
        print('Num. atoms: ' + str(len(self.get_atoms())))
        print('Num heavy atoms: ' + str(self.num_heavy_atoms))
        print('Num. bonds: ' + str(len(self.get_bonds())))
        print('Formula: " ' + self.formula_str + '"')
        print('Mass: ' + str(self.mass))

    def print_atom_bonds(self):
        atoms = self.atoms
        for an_atom in atoms:
            print(an_atom.get_info())
        bonds = self.bonds
        for a_bond in bonds:
            print(a_bond.get_info())


# reads mol v2000
def load_sdf(fpath, ftype):
    mol = molecule()
    if not os.path.exists(fpath):
        print('file does not exist: ' + fpath)
        return mol
    fin = open(fpath, 'r')
    started = 0
    num_seen_atoms = 0
    num_seen_bonds = 0
    for a_line in fin:
        if started == 1:
            if num_seen_bonds < num_bonds and num_seen_atoms == num_atoms:
                num_seen_bonds += 1
                from_index = int(a_line[0:3])
                to_index = int(a_line[3:6])
                bond_type = int(a_line[6:9])
                mol.add_a_bond(from_index, to_index, bond_type)
            if num_seen_atoms < num_atoms:
                num_seen_atoms += 1
                x = float(a_line[0:10])
                y = float(a_line[11:20])
                z = float(a_line[21:30])
                atom_name = a_line[31:35].replace(' ', '')
                if atom_name.lower() == "ac" or atom_name.lower() == "ph" or atom_name.lower() == "me" or atom_name.lower() == "mo":
                    atom_name = "H"
                mol.add_an_atom(atom_name, x, y, z)
        if 'V2000' in a_line:
            num_atoms = int(a_line[0:3])
            num_bonds = int(a_line[3:6])
            started = 1
    fin.close()

    mol.calculate_formula()
    mol.set_mass()
    mol.calculate_networkx()
    mol.calculate_num_protons_attached_to_heavy_atoms()
    mol.set_nghs()
    return mol
