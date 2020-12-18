from collections import defaultdict
import numpy as np

class Logger:
    """Class for logging geometries and energies from a dynamics run.

    Internally, Logger is a dictionary which takes input in the form of a tuple containing
    a string and the thing to be logger.

    Logger is initialized with keys and filenames where the values associated with these
    keys should be printed. Values will be passed in from outside programs.

    After a value is logged, it will be discarded.
    """
    def __init__(self, keys_and_files: list, logging_stride: int, atom_labels=None):
        self.keys, self.files = map(list, zip(*keys_and_files))
        self.data_log = defaultdict(list)
        self.file_log = dict(zip(self.keys, self.files))
        self.file_log.pop('extras', None)

        self.logging_stride = logging_stride
        self.logging_counter = 0

        # atoms for logging the special logging functions
        self.atom_labels = atom_labels

        # check the type of the value and log appropriately
        self.special_logging_functions = {
            np.ndarray: self.log_matrix,
            np.float64: self.log_float
        }

        self.special_logging_keys = [
            "geometry",
            "velocity",
            "force",
        ]

        # initialize all of the non-special files for output
        keys_by_file = defaultdict(list)
        for key, file in self.file_log.items():
            if key not in self.special_logging_keys:
                keys_by_file[file].append(key)
            else:
                # just open the file so it deletes anything that used to be there
                    with open(file, 'w') as f:
                        pass

        
        for file, key_label in keys_by_file.items():
            with open(file, 'w') as f:
                for i, label in enumerate(key_label):
                    f.write(f"Column {i+1}: {label}\n")
                f.write('\n')

    def log(self):
        """Increments the logging counter and possibly calls self.log_data() to
        write the actual data.

        TODO: Create a batch size for the logging so that the actual writing of data occurs in chunks.
        """
        self.logging_counter += 1
        if self.logging_counter >= 10 * self.logging_stride:
            self.log_data()
            self.reset_logging()
    
    def log_data(self):
        """Calls the logging function associated with any special keys.
        Otherwise it just prints it to the appropriate file without any formatting.
        """
        values_by_file = defaultdict(list)
        for key, value in self.file_log.items():
            if key in self.special_logging_keys:
                self.special_logging_functions[type(self.data_log[key][0])](key)
            else:
                values_by_file[value].append(self.data_log[key])
        
        for file, all_data in values_by_file.items():
            zipped_data = list(zip(*all_data))
            with open(file, 'a') as f:
                for values in zipped_data:
                    for value in values:
                        f.write(f" {value:.7f} ")
                    f.write('\n')

    def log_float(self, key):
        """
        Special method for logging floats. This is basically to record the many-body energies, but is also a default logger.
        """
        values = self.data_log[key]
        with open(self.file_log[key], 'a') as f:
            for value in values:
                f.write(f" {value:.7f}\n")


    def log_matrix(self, key):
        """Special method for logging a matrix. Each row is assumed to correspond to the
        atoms in self.atom_labels if they have the same lengths. Otherwise, we just write the
        matrices.

        Args:
            key (string): key to access output file and output data
        """
        values = self.data_log[key]
        if self.atom_labels is not None:
            with open(self.file_log[key], 'a') as f:
                for matrix in values:
                    assert(len(matrix) == len(self.atom_labels))
                    output = str(len(matrix)).rstrip('\n') + '\n' + '\n'
                    for i, value in enumerate(matrix):
                        output += (self.atom_labels[i] + " " +
                                np.array2string(value, precision=14, separator=' ', suppress_small=True).strip('[]') + '\n')
                    f.write(output)
        else:
            with open(self.file_log[key], 'a') as f:
                for value in values:
                    f.write(str(value).strip('[]'))

    def reset_logging(self):
        self.logging_counter = 0
        self.data_log = defaultdict(list)