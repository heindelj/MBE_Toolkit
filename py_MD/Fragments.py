import numpy as np
from tempfile import NamedTemporaryFile
import itertools

class Fragments:
    def __init__(self, xyz_file):
        self.xyz_file = xyz_file
        self.header, self.atom_labels, self.fragments = self.get_fragments_from_xyz_file()
    
    def fragment_geometry(self, geometry):
        """Takes an array of cartesian coordinates and splits it into fragments
        according to the shape of self.fragments.

        Args:
            geometry (ndarray): Nx3 array of cartesian coordinates
        """
        self.fragments = np.reshape(geometry, self.fragments.shape)

    def write_fragment_to_temporary_file(self, fragment, array_indices):
        """Takes a fragment which is just the matrix of xyz coordinates and prints them 
        out to a file which can be read by TTM2.1-F. Also takes the array indices
        to identify the appropriate atom labels.
        """

        output = str(len(fragment)) + '\n\n'
        for i in range(len(fragment)):
            output += (self.atom_labels[array_indices[i]] + " " +
                    np.array2string(fragment[i], precision=14, separator=' ', suppress_small=True).strip('[]') + '\n')
 
        temp_file = NamedTemporaryFile('w', delete=False)
        with open(temp_file.name, 'w') as f:
            f.write(output)
        return temp_file

    def get_indices_for_fragment_combination(self, i_order: int):
        """Takes the order of MBE we're doing and returns a list of tuples containing
        the indices into the original array of fragments. Also returns a list of tuples
        containing the atom indices for each fragment, so that we can index into both
        the fragments and atoms of the total system.

        This only has to be done once as long as the potential can guarantee to return
        the forces in the same order as atoms are given to the potential.

        Args:
            i_order       (int): order of the mbe we're currently working on
        """
        fragment_index_array = [x for x in range(0, len(self.fragments))] # e.g. [0, 1, 2, 3]
        combinations = list(itertools.combinations(fragment_index_array, i_order))

        # this tells you the number of atoms which appear before the nth fragment begins,
        # so basically the size of all preceding fragments. This is kind of ugly, but
        # I don't see a simpler way to do this.
        distance_to_nth_fragment = []
        distance_until_now = 0
        for i, frag in enumerate(self.fragments):
            if i-1 > -1:
                distance_until_now = distance_to_nth_fragment[i-1]
            distance_to_nth_fragment.append(frag.shape[0] + distance_until_now)
            distance_to_nth_fragment[i] -= distance_to_nth_fragment[0]

        atom_index_array = []
        for fragment_indices in combinations:
            list_of_lists = [list(range(distance_to_nth_fragment[x], distance_to_nth_fragment[x] + len(self.fragments[x]))) for x in fragment_indices]
            list_of_lists = [item for sublist in list_of_lists for item in sublist]
            atom_index_array.append(list_of_lists)

        return atom_index_array

    def make_nmers(self, mbe_order):
        """Returns a list of numpy arrays of all n-mers of order mbe_order.

        e.g. If mbe_order=2, returns a list of all dimers made from self.fragments.fragments

        Args:
            mbe_order (int): Order of the mbe to form nmers of (monomers, dimers, etc.)
        """
        combinations = list(itertools.combinations(self.fragments, mbe_order))
        for i, arrays in enumerate(combinations):
            combinations[i] = np.vstack(tuple(arrays))
        return combinations

    def get_fragments_from_xyz_file(self):
        """Reads an xyz file containing a single geometry where fragments are delimited by '--'.

        Input: string representing path to input file
        Returns: header of xyz file, atom labels, and numpy array of numpy arrays of xyz coordinates of each fragment
        """
        fragments = []
        atomLabels = []
        header = []
        with open(self.xyz_file) as ifile:
                line = ifile.readline().split()
                if line and line[0].isdigit():
                    natoms = int(line[0])
                    title = ifile.readline()
                    try:
                        header.append(str(natoms) + '\n' + ' '.join(line[1:]) + title)
                    except IndexError:
                        header.append(str(natoms) + '\n' + title)

                    while True:
                        fragment__ = []
                        # get first line of coordinates
                        coord_line = ifile.readline()
                        if not coord_line:
                            break
                        while '--' not in coord_line and coord_line:
                            line = coord_line.split()
                            atomLabels.append(line[0])
                            fragment__.append(list(map(float, line[1:4])))
                            coord_line = ifile.readline()
                        fragments.append(np.array(fragment__))
        return header, atomLabels, np.array(fragments, dtype=np.float64)

    def write_geoms(self, optional_output="", ofile=None):
        """
        args: header, labels, and coords as output by read_geoms
        return: no return

        Writes the the molecules to stdout in xyz format if no ofile is specified.
        Otherwise, write the geometries to ofile.
        """
        output = str(self.header[0]).rstrip('\n') + '\n' + str(optional_output) + '\n'
        for i in range(len(self.fragments)):
            for j in range(len(self.fragments[i])):
                output += (self.atom_labels[i*len(self.fragments[i]) + j] + " " +
                        np.array2string(np.array(self.fragments[i][j]), precision=14, separator=' ', suppress_small=True).strip('[]') + '\n')
        if ofile is None:
            print(output)
        else:
            with open(ofile, 'w') as f:
                f.write(output)
    
    def write_single_geometry(self, geometry, ofile=None):
        output = str(len(geometry)).rstrip('\n') + '\n' + '\n'
        for i in range(len(geometry)):
            output += (self.atom_labels[i] + " " +
                    np.array2string(np.array(geometry[i]), precision=14, separator=' ', suppress_small=True).strip('[]') + '\n')
        if ofile is None:
            print(output)
        else:
            with open(ofile, 'w') as f:
                f.write(output)