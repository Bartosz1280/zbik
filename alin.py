# External imports
import os
import re
import itertools
import numpy as np
from typing import *
from datetime import datetime
# Internal imports
from loader import *
from sequence import *

class Alin:
    """
    Class holding alignments
    """
    def __init__(self):
        # Script related
        self.path = None
        self.report = None # Holds strings for making a text file with a report
        self.seq_type = None # Holds info about sequences type within the instance
        self.head_msg = str(f">> zbik.Alin {datetime.now().strftime('%H:%M')} :")

        # DNA related
        if self.seq_type == "Nucleotides": # Initiated DNA/RNA specific attributes
                self.point_mutations_positions = None # list of mutated nucleotides position
                self.point_mutations_number = None # TBI : calculates number of points mutation
                self.undetermined_mutations_number = None # TBI : number of cases where WT couldn't be determined
                self.__deter__ = None # Secret attribute: see: determine_point_mutations_positions()
                self.tt_ratio = None # Transitions and Transversions ratio
                # To be implemented
                self.tt_positions = None # Transitions and Transversions positions

        # Common for sequences
        if self.path:
            self.sequences = loader(self.path)  # Automatically cast sequences when path provided
        else:
            self.sequences = None
        self.__seqs_type_conflict = None # Invokes when more than one sequence type was loaded
        self.distance_matrix = None  # Holds generated  distance matrix
        self.edit_distances_matrix = None
        # Attributes below were not implemented yet
        self.matrix = None # Holds substitution matrix
        self.scores = None
        self.master = None

        print(f"{self.head_msg} Alin instance initiated")
    def __str__(self):
        sum_str= ">>> Alin instance >>>\n" \
                 " Bounded Sequences:\n"
        for seq in enumerate(self.sequences):
            ind, seq_obj = seq
            if ind < 10:  # Nicer index placement for ind > 10
                sum_str=sum_str + f" - {ind} : {seq_obj.sequence_type}: {seq_obj.instance_name}\n"
            else:
                sum_str=sum_str + f" - {ind}: {seq_obj.sequence_type}: {seq_obj.instance_name}\n"
        return sum_str
    def show_point_mutations_positions(self):
        def head_mark_positions(indexes,Seq):
            seq_str = list(Seq.sequence)
            for i in indexes:
                seq_str[i] = '\x1b[1;37;41m'+seq_str[i]+'\x1b[0m'
            return Seq.instance_name + ' : ' +"".join(seq_str)
        # Master function execution
        print(f"{self.head_msg} Position of point mutations")
        for seq in self.sequences:
                print(head_mark_positions(self.point_mutations_positions,seq))

    # Functions associated with building up and use of the instance
    def bounded_sequences(self):
        """
        Class Alin method.

        Function for printing out Sequence objects bounded to Alin in user-friendly mode
        [index] [accession number] [sequence type]
        """
        for seq in enumerate(self.sequences):
            ind, seq_obj = seq
            if ind < 10: #Nicer index placement for ind > 10
                print(f" - {ind} : {seq_obj.sequence_type}: {seq_obj.instance_name}")
            else:
                print(f" - {ind}: {seq_obj.sequence_type}: {seq_obj.instance_name}")
    def help(self):
        ...
    def display_matrix(self,matrix):
        match matrix.lower():
            case 'distance' | 'distance_matrix':
                matrix = self.distance_matrix
            case 'edit_distances' | 'edit_distances_matrix':
                matrix = self.edit_distances_matrix
        for r in matrix:
            print("".join("".join(str(r).split("[")).split("]")) )
    def cast(self, path=None, seqs=None):
        """
        WHEN path is provided:
            Class Alin build-in function. Load sequences from files into the instance.
            Created Sequence instances are in globals, whereas the Alin instance links to their attributes
        WHEN path is not provided:
            Takes already created Sequence class instances in seqs argument.
        :param seqs: instance(s) of Sequence class
        :param path: str() of path to directory with sequences files in .fasta/.gb/.gp formats
        """

        def loader(dir_path):
            """
            Loads multiple DNA/RNA/protein sequences files in .fasta, .gb ot .gp formats into separate Sequence class instances.
            :param dir_path: path to directory where files are stored
            """
            inited_seq_ids = list()  # Holds ids of Sequence instances created by the function
            #  Messages used later for prints
            invalid_file_msg = 'Invalid file'
            head_msg = str(f">> zbik.loader {datetime.now().strftime('%H:%M')} :")

            def file_loader(path):
                """
                Loads a NCBI file. Detects if the file is in .fasta or .gb format
                Detect if the pass sequence is a nucleotide or amino acis sequence
                :param path: (str) of path to NCBI file
                :return:
                """
                def sequence_filter(file_slice):
                    """
                    Sub-function of file_loader(). Holds an expression which enables to get a sequence from a long string,
                    produced by splicing the input text file
                    :param file_slice: str from a list produced with file slicing (i.e. Regex splice on ORIGIN in .gb files)
                    :return seq: str of a sequence
                    """
                    seq = "".join(
                        list(filter(lambda x: x in ['a', 'c', 't', 'g', 'n', 'u', 'h', 'v', 'b', 'r', 'd', 'f', 'e',
                                                    'z', '+', 'x', 'j', 'm', 'i', 'q', 'p', 'y', 'k', 'l', 's',
                                                    'w', '-'], file_slice.lower()))).upper()
                    return seq

                with open(path) as file:
                    file = file.readlines()
                # Detecting file type and get the sequence
                match file[0][0]:
                    case ">":  # .fasta file case
                        file_type = 'FASTA'
                        acs_number = file[0].split(" ")[0][1::]  # Splits the first line to get accession number
                        seq = sequence_filter("".join(file[1:]))  # Takes everything from the second line to the last
                        # and filters the file slice (str) using sequence_filter() function
                    case "L":  # .gb file case
                        if file[0][0:5] == "LOCUS":  # Detects LOCUS - the first word in the .gb file
                            file_type = 'GeneBank'
                            for i in range(0, len(file)):
                                # Gets the accession number with version sequences string from .gb/ .gp file
                                if list(filter(lambda x: x != '', file[i].split("\n")[0].split(" ")))[0] == "VERSION":
                                    acs_number = list(filter(lambda x: x != '', file[i].split("\n")[0].split(" ")))[1];
                                    break
                                if i == len(file) - 1:
                                    print("Not able to determine accession number from the provided GeneBank file")
                            seq = sequence_filter(
                                re.split("ORIGIN", "".join(file))[1])  # Takes everything form index 1 after
                            # slicing a file on ORIGIN string and filters the file slice (str) using sequence_filter() function
                        else:
                            print(invalid_file_msg);
                            return
                    case _:
                        print(invalid_file_msg);
                        return
                return acs_number, file, file_type, seq  # Tuple used by the next function and unpacked there

            def instance_builder(file_loader):
                """
                Initiate a Sequence class instance
                :param file_loader: (tuple) output of file_loader() function
                """
                acs_number, file, file_type, sequence = file_loader  # Unpacks the output tuple from file_loader
                instance_name = acs_number  # Define name of a new instance that will be created using accession number
                if "." in instance_name: instance_name = "_".join(
                    instance_name.split("."))  # Removes dots from acs_number
                instance_count = 0  # counting for creating alternative instance name
                while instance_name in globals():  # Starts when the instance already exists
                    print(f"{head_msg} Existing {instance_name} instance name. Creating a new name")
                    instance_count += 1
                    instance_name = instance_name + "_v" + str(instance_count)  # Creates alternative instance name
                # Create a new sequence instance
                globals()[instance_name] = Sequence(acs_number, sequence, file, file_type, instance_name)
                print(f"{head_msg} {instance_name} Sequence instance initiated")
                inited_seq_ids.append(instance_name)

            # Loop for loading multiple sequence files from a given directory
            for path in list(map(lambda x: dir_path + "/" + x, os.listdir(dir_path))):
                instance_builder(file_loader(path))
            return inited_seq_ids
        ### Master function execution
        if path and seqs:  # Handle sequences source conflicts
            print(f"{self.head_msg} You provided two possible sources of sequences files.")
            print(">>> Which one you want to use?\n (1) - directory path \n (2) - Sequence instances")
            inp = input("")
            if inp == 1:
                seqs = None
                print(">>> Sequences purged purged from the instance")
            elif inp == 2:
                path = None
                print(">>> Path purged from the instance")
        if not path and seqs:  # Use pre-existing instances
            if len(seqs) > 1:
                self.sequences = seqs
                print(f"{self.head_msg} Number of added sequences: {len(seqs)} ")
            else:
                print(f"{self.head_msg} List of sequences is empty!")
                return
        if path and not seqs:  # Loads Sequence instances from files
            self.sequences, cnt = list(), 0
            seq = loader(path)
            for seq_id in seq:
                cnt += 1
                self.sequences.append(globals()[seq_id])
            print(f"{self.head_msg} Number of added sequences: {cnt} ")
            self.__seqs_type_conflict = self.__detect_seq_type_conflict()
    def __detect_seq_type_conflict(self):
        """
        :return: bool() True if there are more than 1 sequences types
        """
        try:
            if 1 == len(set(list(map(lambda x: x.sequence_type, self.sequences)))):
                return False
            else:
                print(f"{self.head_msg} Conflicting sequences types detected!")
                return True
        except TypeError:
            if not seq:
                print(f"{self.head_msg} Cannot check for conflicting sequences types. Instance is not initialized!")
                return None
            else:
                if len(self.seq) == 0:
                    print(f"{self.head_msg} : Cannot check for conflicting sequences types. Sequence list is empty!")
                    return None
    # Nucleotide type specific methods
    def generate_edit_distances_matrix(self):
        def wagner_fischer(str1, str2):
            """
            Implementation of simplified Wagner-Fischer algorithm
            :param str1: First str
            :param str2: Second str
            :return: int of operations needed to transform one string to another
            """
            d = np.array(np.zeros((len(str1) + 1) * (len(str2) + 1))).reshape((len(str1) + 1), (len(str2) + 1))
            for r in range (d.shape[0]): # Number of rows
                for c in range(d.shape[1]): # Number of columns
                    if r == 0: # row with 0 index
                        d[r][c] = c
                    elif c == 0: # column with 0 index
                        d[r][c] = r
                    elif str1[r - 1] is str2[c - 1]: # letter on index of current column -1 == letter on index of current row -1
                        d[r][c] = d[r - 1][c - 1]
                    else:
                        d[r][c] = 1 + min(d[r][c - 1], d[r - 1][c], d[r - 1][c - 1])
            return int(d[len(str1)][len(str2)])
        edit_distances=[]
        for combo in itertools.product(self.sequences, repeat=2):
            edit_distances.append(wagner_fischer(combo[0].sequence, combo[1].sequence))
        self.edit_distances_matrix = np.matrix(edit_distances).reshape(len(self.sequences), len(self.sequences))
        print(f"{self.head_msg} Matrix of edited distances for {len(self.sequences)} sequences created\n")
        self.display_matrix('edit_distances_matrix')
        print("")
        self.bounded_sequences()
    def determine_point_mutations_positions(self):
        """
        Class Alin method.
        Determines point mutation by assessing bounded sequences by column (position).
        When at least one sequence differs index (position in sequence) is added to a list of point mutation
        position stored in attribute point_mutations_positions .
        :return:
        """

        def detect_mutation_on_nucleotide(aln_tuple):
            """
            Sub function of determine_point_mutations_positions().
            Detects mutations in a column (nucleotide position).
            :param aln_tuple: column of nucleotides or AA
            :return: True when more than one possible nucleotide presents in the column
            """
            if len(set(aln_tuple)) != 1:
                return True
            else:
                return False

        def determine_point_mutations_quantity(self):
            """
            Sub function of determine_point_mutations_positions().
            Determines possible point mutations quantity within bounded sequences.
            :param self:
            :return:
            """
            def count_possible_mutations(aln_column):

                nucleo_dict = dict().fromkeys(aln_column)
                for key in nucleo_dict:
                    nucleo_dict[key] = list(aln_column).count(key)
                # Conditional statement below check if the highest value occurs only once in the dictionary
                if list(nucleo_dict.values()  # Creates the list of nucleo_dict values
                        ).count(nucleo_dict[max(nucleo_dict)]# Counts how many times the highest value occurs in the lust
                                ) > 1:  # Checks if the highest value in nucleo_dict occurs only once
                    return False  # Returns False when cannot determine WT. False acts like 0, but can be captured with conditional statement
                else:  # Condition for only one possible WT phenotype
                    if len(nucleo_dict.keys()) > 1:  # Conditional that's prevents returning of negative int's
                        return len(nucleo_dict.keys()) - 1  # Number of possible mutation (all possible nucleotides - WT)
                    else:
                        return 0  # Returns 0 if subtraction would give -1

            # Execution of the master function
            self.__deter__ = list( # Creates a booleans list with detect_mutation_on_nucleotide()
                map(count_possible_mutations, np.transpose( # Transpose alignment array to access columns
                    list(map(lambda x: list(x.sequence), self.sequences # Converts string of sequences to list of letters
                             )))))
            # Assignment below count number of cases when it was not possible to determine wild type nucleotide/AA
            self.undetermined_mutations_number = sum(list(filter( #Gets number of False occurrence
                lambda x: x == False, self.__deter__))) # Filters False values from a list of False and int of possible mutations
            self.point_mutations_number = sum(self.__deter__) # Sums number of possible mutations: False = 0 in Python

        # Execution of the master function
        self.point_mutations_positions=list( # return list of indexes to point_mutation attribute
            map(lambda x : x[0], # Takes only indexes of point mutation using integer from tuple
            filter(lambda x : x[1] == True, # Filters only detected point mutations using boolean from tuple
            enumerate(map( # Enumerates to gets index, boolean (point mutation) tuples list
            detect_mutation_on_nucleotide, np.transpose( # Checks for point mutations on each column in alignment
                    list(map(lambda x: list(x.sequence), self.sequences #  Created an alignment by building list of
                             # nucleotides in each sequence, assembling it into np.array, transposing it so the
                             # map/loops access columns
                             ))))))))

        print(f"{self.head_msg}  point mutations detected on {len(self.point_mutations_positions)} different positions")
        determine_point_mutations_quantity(self)

    def get_tt_ratio(self):
        """
        Class ALin method.
        Calculates Transition/Transversion ratio
        :return:

        """
        mut = [
            np.transpose(list(map(lambda x: list(x.sequence), self.sequences)))[ind] for ind in self.point_mutations_positions
        ]
        if len(self.sequences) > 2:
            print('More than 2 sequences will be implemented soon')
            return None
        else:
            def determine_transition(aln_column):
                """
                :param aln_column:
                :return:
                """
                col = [aln_column[0],aln_column[1]] # List declaration to prevent alternation in original sequences
                col.sort()                          # with sort() function
                col = "".join(col)
                match col:
                    case "AG": # transition A<->G
                        return True
                    case "CT": # transition A<->G
                        return True
                    case _: # transversions
                        return False
        tt_list = list(map(determine_transition,mut))
        self.tt_ratio=tt_list.count(True)/tt_list.count(False)
        print(f"{self.head_msg} transition/transversion ratio R equals {round(self.tt_ratio,3)}")

    def generate_distance_matrix(self):
        distances= list()
        for combo in itertools.product(self.sequences,repeat =2):
            distances.append(len(list(
                filter(lambda x: len(x) != 1, map(
                    lambda x: set(x), np.array(
                        [list(combo[0].sequence), list(combo[1].sequence)]).transpose())))) / len(
                combo[0].sequence))

        self.distance_matrix = np.matrix(distances).reshape(len(self.sequences),len(self.sequences))
        print(f"{self.head_msg} Distance matrix for {len(self.sequences)} sequences created\n")
        self.display_matrix("distance")
        print("")
        self.bounded_sequences()
