from datetime import datetime
from sequence import *
import re
import os


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
            seq = "".join(list(filter(lambda x: x in ['a', 'c', 't', 'g', 'n', 'u', 'h', 'v', 'b', 'r', 'd', 'f', 'e',
                                                      'z', '+', 'x', 'j', 'm', 'i', 'q', 'p', 'y', 'k', 'l', 's',
                                                      'w', '-'], file_slice.lower()))).upper()
            return seq
        with open(path) as file: file = file.readlines()
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
                    for i in range(0,len(file)):
                        # Gets the accession number with version sequences string from .gb/ .gp file
                        if list(filter(lambda x : x != '',file[i].split("\n")[0].split(" ")))[0] == "VERSION":
                            acs_number = list(filter(lambda x : x != '',file[i].split("\n")[0].split(" ")))[1]; break
                        if i == len(file)-1:
                            print("Not able to determine accession number from the provided GeneBank file")
                    seq = sequence_filter(re.split("ORIGIN", "".join(file))[1])  # Takes everything form index 1 after
                    # slicing a file on ORIGIN string and filters the file slice (str) using sequence_filter() function
                else: print(invalid_file_msg); return
            case _: print(invalid_file_msg); return
        return acs_number, file, file_type, seq  # Tuple used by the next function and unpacked there

    def instance_builder(file_loader):
        """
        Initiate a Sequence class instance
        :param file_loader: (tuple) output of file_loader() function
        """
        acs_number, file,  file_type, sequence = file_loader  # Unpacks the output tuple from file_loader
        instance_name = acs_number  # Define name of a new instance that will be created using accession number
        if "." in instance_name: instance_name = "_".join(instance_name.split("."))  # Removes dots from acs_number
        instance_count = 0  # counting for creating alternative instance name
        while instance_name in globals():  # Starts when the instance already exists
            print(f"{head_msg} Existing {instance_name} instance name. Creating a new name")
            instance_count += 1
            instance_name = instance_name+"_v"+str(instance_count)  # Creates alternative instance name
        # Create a new sequence instance
        globals()[instance_name] = Sequence(acs_number, sequence, file, file_type, instance_name)
        print(f"{head_msg} {instance_name} Sequence instance initiated")
        inited_seq_ids.append(instance_name)
    # Loop for loading multiple sequence files from a given directory
    for path in list(map(lambda x: dir_path + "/" + x, os.listdir(dir_path))):
        instance_builder(file_loader(path))
    return inited_seq_ids
