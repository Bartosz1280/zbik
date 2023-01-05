import sys
import platform
from alin import *
import numpy as np
run = True

green_tick = f"\33[32m✓\033[0m"
red_cross = f"\033[91m✕\033[0m"
def display_file_message(message):
    path = f"messages/{message}"
    with open(path) as msg:
        msg = msg.read()
    print("")
    print(msg)
    print("")
def print_bold(msg):
    print("\033[1m" + msg + "\033[0m")
display_file_message("logo")
# Checks for compatibility
if int(sys.version.split(" ")[0].split(".")[1]) <10:
    print(f"{red_cross} : Python version  is {sys.version.split(' ')[0]}")
    display_file_message("python_warn.txt")
    run = False
else:
    print(f"{green_tick} : Sufficient Python version")
if platform.platform().split("-")[0] != 'Linux':
    display_file_message("os_warn.txt")
else:
    print(f"{green_tick} : OS compatible with one used for testing.")
if int(np.__version__.split(".")[1]) < 23:
    print(red_cross+" : You are using older version of NumPy  (<1.23)")
else:
    print(green_tick+" : Sufficient NumPy version")
    print_bold("\nWelcome to Zbik 0.0.1. Program is in early development phase")
    print_bold("Some features will be polished soon...")
while run:
    print_bold("\n>What would you like to do:")
    print("  1 - Read about current version")
    print("  2 - See how to use Zbik classes")
    print("  3 - Run Rosalind examples")
    print("  E - exit")
    match input("? ").upper():
        case "E":
            print("See you next time!")
            run = False
        case "1":
            display_file_message("current_release.txt")
            input("Press enter to continue ...")
        case '2':
            man = True
            while man:
                print_bold("\n>Which one would you like to see?:")
                print("  1 - Class Sequence")
                print("  2 - Class Alin")
                print("  3 - Other functionality")
                print("  E - return to main menu")
                match input("? ").upper():
                    case "E":
                        man = False
                    case "1":
                        print("\nBy providing a path to function loader(), instances of Sequences")
                        print("are created with names associated to their accession numbers\n")
                        print(">> loader('./sequences')")
                        input("Press enter to continue ...\n")
                        seq_instance = loader("./sequences")
                        input("Press enter to continue ...\n")
                        print("\nIn case you want to check initialized Sequences you can use")
                        print("ids attribute of any loaded Sequence instance. Also loader function\n"
                              "returns a list of instances names. This method of loading sequences\n"
                              "is not recommended. Loading with class Alin instance is preferred.")
                        print("\nTo get a quick access to basic information about the instance")
                        print("You can just call a print() function on the instance\n")
                        input("Press enter to continue ...\n")
                    case "2":
                        print("\nClass Alin holds alignment of sequences. For now it is a main\n"
                              "class of Zbik, that holds most of its functionality. You can call it\n"
                              "without any attributes. You will get an information that the Alin \n"
                              "instance was initialized\n")
                        print(">>> Test_Alin = Alin()")
                        input("Press enter to continue ...\n")
                        Test_Alin = Alin()
                        input("Press enter to continue ...\n")
                        print("\nSince the Alin instance is empty, all its attributes are assigned as None\n"
                              "You can load multiple sequences from one directory using cast() method.\n")
                        print(">>> Test_Alin.cast('./sequences')")
                        Test_Alin.cast('./sequences')
                        print("\nYou might notice that Alin instance warned you about conflicting sequences.\n"
                              "This occurs when different types of sequences are loaded into the same Alin instance\n"
                              "Lets examine what sequences are loaded into Alin instance. It can be done using \n"
                              "bounded_sequences() method of the class\n")
                        input("Press enter to continue ...\n")
                        Test_Alin.bounded_sequences()
                        print("\nAs you can see Test_Alin holds two DNA and two amino acids sequences\n")
                        input("Press enter to continue ...\n")
                        print("You can also call a print function on Alin instance to get \n"
                              "information about sequences it holds. In the future this methods\n"
                              "will display much more details about the instance")
                        print(">>> print(Test_Alin)")
                        input("Press enter to continue ...\n")
                        print(Test_Alin)
                        input("Press enter to continue ...\n")
                    case "3":
                        print("\nTo see examples of other functionality check:")
                        print(" Run Rosalind examples in main menu")
                        input("Press enter to continue ...\n")
        case '3':
            rosalind = True
            while rosalind:
                print_bold("\n>Which one would you like to see?:")
                print("  1 - Counting Point Mutations")
                print("  2 - Transitions and Transversions")
                print("  3 - Creating a Distance Matrix")
                print("  4 - Creating a edit distances")
                print("  E - return to main menu")
                match input("?").lower():
                    case "e":
                        rosalind = False
                    case "1":
                        display_file_message("rosalin_counting_point_mutations.txt")
                        input("Press enter to continue ...\n")
                        print("")
                        print("Alin class have its method called: determine_point_mutations_positions()\n"
                        "Another method: show_point_mutations_positions() shows positions \n"
                        "Of point mutations in human readable version.\n"
                        "For demonstration purpose the example data set of the problem was used")
                        print(">>> Point_Mutations_Alin.determine_point_mutations_positions()")
                        input("Press enter to continue ...\n")
                        Point_Mutations_Alin = Alin()
                        path = './rosalind/point_mutations'
                        Point_Mutations_Alin.cast(path)
                        print("")
                        Point_Mutations_Alin.determine_point_mutations_positions()
                        print("\nPoint mutations were determined on 7 different positions.")
                        print("The instance stores them in form of list of indexes")
                        print("Using show_point_mutations_positions(). Shows it in human readable form\n")
                        print(">>> Point_Mutations_Alin.show_point_mutations_positions()")
                        input("Press enter to continue ...\n")
                        Point_Mutations_Alin.show_point_mutations_positions()
                        input("Press enter to continue ...\n")
                    case "2":
                        display_file_message("rosalin_tt_ratio.txt")
                        input("Press enter to continue ...\n")
                        print("\n Alin class has a method called get_tt_ratio(). Results of the")
                        print("function are stored in tt_ratio attribute\n")
                        print("For demonstration purpose the example data set of the problem was used")
                        input("Press enter to continue ...")
                        TT_alin = Alin()
                        path = './rosalind/transition_transversions'
                        TT_alin.cast(path)
                        print("\nTo get tt_ratio first you need to run determine_point_mutations_positions()")
                        print("method. To speed up processing Alin will ignore conserved nucleotides during")
                        print("iteration for get_tt_ratio() \n")
                        input("Press enter to continue ...\n")
                        TT_alin.determine_point_mutations_positions()
                        TT_alin.get_tt_ratio()
                        print(f"Script determined R equals to {TT_alin.tt_ratio:.3f}")
                        input("Press enter to continue ...\n")
                    case "3":
                        display_file_message('rosalin_distance_matrix.txt')
                        input("Press enter to continue ...\n")
                        print("\nAlin class has its method called generate_distance_matrix()")
                        print("It's allows to obtain a distance matrix of all loaded sequences")
                        print("Results are stored in attribute distance_matrix attribute.")
                        print("When the calculation are completed Alin prints out the results.\n")
                        print("For demonstration purpose the example data set of the problem was used")
                        input("Press enter to continue ...\n")
                        DM_Alin = Alin()
                        path='./rosalind/distance_matrix'
                        DM_Alin.cast(path)
                        DM_Alin.generate_distance_matrix()
                        print("\n NOTE: Obtained matrix might look a bit different from the sample output")
                        print("This is cause by the order at which Alin loads files\n")
                        input("Press enter to continue ...\n")
                    case "4":
                        path='./rosalind/rosalin_edit_distance.txt'
                        display_file_message('rosalin_edit_distance.txt')
                        input("Press enter to continue ...\n")
                        print("Alin class has method called generate_edit_distances_matrix()")
                        print("This function determines edited distances between all passed sequences")
                        print("Result is presented as a matrix and its stores in edit_distances_matrix")
                        input("Press enter to continue ...\n")
                        ED_alin = Alin()
                        path ='./rosalind/edit_distances'
                        ED_alin.cast(path)
                        ED_alin.generate_edit_distances_matrix()
                        print("\n NOTE: Returning edited distances as a matrix allows to check more than two sequences")
                        input("Press enter to continue ...\n")