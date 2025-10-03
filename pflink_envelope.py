import os
import sys
import argparse
from rich_argparse import RichHelpFormatter

import argparse
import os
import sys

from pigeon_feather.hxio import *
from rich_argparse import RichHelpFormatter

from HXMS_IO import write_hxms_file
from Helper_Functions import _conver_PFhxms_to_hxms
from Parsers import FlagsParser, _parse_dynamX, _parse_HDXWorkbench, _parse_biopharma, _parse_byos, _parse_HDExaminer, \
    _parse_custom

# ANSI escape codes for colors
RED = "\033[91m"
RESET = "\033[0m"

def main():
    # Set up the command-line argument parser with the custom formatter


    parser = argparse.ArgumentParser(
        description="Process HX-MS data.",
        formatter_class=lambda prog: RichHelpFormatter(prog, max_help_position=200, width=400)
    )

    parser = argparse.ArgumentParser(
        description="Process raw H/D exchange mass spectrometry data into HXMS file format.")


    #I/O and File Arguments
    parser.add_argument("--input_csv_path", type=str, required=True,
                        help="Path to the file containing the raw input CSV/data files.")
    parser.add_argument("--output_hxms_path", type=str, required=True,
                        help="Path to the output directory where the final .hxms files will be saved.")
    parser.add_argument("--flags_file_path", type=str, required=False,
                        help="Path to a flags or configuration file to override default settings.")
    parser.add_argument("--peptide_list", type=str, required=False,
                        help="Path to a CSV/TXT file containing a list of peptides to filter or process.")
    parser.add_argument("--raw_spectra_path", type=str, required=True,
                        help="Name of the HDExam spectra folder inside the input path.")

    #Metadata Arguments
    parser.add_argument("--saturation", type=float, required=False,
                        help="The $\text{D}_2\text{O}$ saturation percentage (0.0 to 1.0) to record in the file header.")
    parser.add_argument("--ph", type=float, required=False,
                        help="The measured $\text{pH}$ ($\text{pD}$ read) value of the exchange buffer.")
    parser.add_argument("--temperature", type=float, required=False,
                        help="The temperature of the exchange reaction in Kelvin ($\text{K}$).")
    parser.add_argument("--protein_name", type=str, required=False,
                        help="The common or unique name of the protein.")
    parser.add_argument("--protein_state", type=str, required=False,
                        help="The experimental state of the protein (e.g., 'Apo', 'Bound', 'Mutant').")
    parser.add_argument("--protein_sequence", type=str, required=False,
                        help="The full amino acid sequence of the protein.")

    #Filter/Choice Argument
    parser.add_argument("--include_exclude",
                        type=str,
                        required=False,
                        choices=['include', 'exclude'],  # Enforces a specific value
                        default=None,
                        help="Determines if the peptide list file should be used to include or exclude peptides from processing.")
    parser.add_argument("--file_type",
                        type=str,
                        required=False,
                        choices=["HDExaminer"],  # Enforces a specific value
                        default=None,
                        help="Input CSV format type")

    # Parse the arguments
    args = parser.parse_args()

    # Construct the full paths from the parsed arguments
    INPUT_PATH = args.input_csv_path
    OUTPUT_PATH = args.output_hxms_path
    FLAGS_PATH = args.flags_file_path
    PEPTIDE_LIST_PATH = args.peptide_list
    RAW_SPEC_PATH = args.raw_spectra_path

    # Path Validations for input files
    if not os.path.isfile(INPUT_PATH):
        print(f"{RED}Error: Your input path '{INPUT_PATH}' does not exist or is not a file! Please correct it and run again.{RESET}")
        sys.exit(1)

    if not os.path.isdir(OUTPUT_PATH):
        print(f"{RED}Error: Your output path '{OUTPUT_PATH}' does not exist or is not a directory! Please correct it and run again.{RESET}")
        sys.exit(1)

    if not os.path.isfile(RAW_SPEC_PATH):
        print(f"{RED}Error: Your spectra path '{RAW_SPEC_PATH}' does not exist or is not a directory! Please correct it and run again.{RESET}")
        sys.exit(1)

    if FLAGS_PATH:
        if not os.path.isfile(FLAGS_PATH):
            print(f"{RED}Error: The table file '{FLAGS_PATH}' does not exist! Please correct it and run again.{RESET}")
            sys.exit(1)
    if PEPTIDE_LIST_PATH:
        if not os.path.isfile(PEPTIDE_LIST_PATH):
            print(f"{RED}Error: The peptide path '{PEPTIDE_LIST_PATH}' does not exist or is not a file! Please correct it and run again.{RESET}")
            sys.exit(1)

    if FLAGS_PATH:
        parser = FlagsParser(FLAGS_PATH)
        flags = parser.parse()

        protein_name_state_info = flags['protein_name_states']
        if protein_name_state_info is None:
            print("You must have a protein name and state in your flags. Example: [['Protein1' 'State1'],['Protein1' 'State2'],['Protein2' 'State3']]")
            exit()
        try:
            if len(protein_name_state_info[0]) != 2:
                print("You must have a protein name and state in your flags. Example: ['Protein' 'State']")
                exit()
        except:
            print("You must have a protein name and state in your flags. Example: ['Protein' 'State']")
            exit()
        protein_name = flags['protein_name_states'][0][0]
        protein_sequence = flags['protein_sequence']
        if protein_sequence is None:
            print("You must have a protein sequence in your flags. Example: ['AAGWDGA']")
            exit()

        saturation = flags['d20_saturation']
        if saturation is None:
            print("You must have a d20_saturation in your flags. Example: 0.6")
            exit()
        try:
            saturation = float(saturation)
            if saturation > 1 or saturation < 0:
                print("Your d20_saturation must be a float between 0 and 1 inclusive. Example: 0.6")
                exit()
        except:
            print("Your d20_saturation must be a float. Example: 0.6")
            exit()

        temperature = flags['temp']
        if temperature is None:
            print("You must have a temp (in K) in your flags. Example: 293.5")
            exit()
        try:
            temperature = float(temperature)
        except:
            print("You must have a temp (in K) in your flags that is an int or float. Example: 293.5")
            exit()
        ph = flags['ph']
        if ph is None:
            print("You must have a ph in your flags. Example: 7.2")
            exit()
        try:
            ph = float(ph)
        except:
            print("You must have a ph in your flags that is an int or float. Example: 7.2")
            exit()
        file_type = flags['file_type']
        if len(flags['protein_sequence']) == 1:
            flags['protein_sequence'] = flags['protein_sequence'][0]
        if file_type is None:
            print("You must choose a file type!")
            exit()


    else:
        saturation = args.saturation
        ph = args.ph
        temperature = args.temperature
        protein_name = args.protein_name
        protein_state = args.protein_state
        protein_sequence = args.protein_sequence
        file_type = args.file_type
        if saturation is None:
            print("You must provide a saturation value!")
            exit()
        if ph is None:
            print("You must provide a ph value!")
            exit()
        if temperature is None:
            print("You must provide a temperature value!")
            exit()
        if protein_name is None:
            print("You must provide a protein_name value!")
            exit()
        if protein_state is None:
            print("You must provide a protein_state value!")
            exit()
        if protein_sequence is None:
            print("You must provide a protein_sequence value!")
            exit()
        if file_type is None:
            print("You must provide a file_type value!")
            exit()
        flags = {
            'protein_name_states': [protein_name + " " + protein_state],
            'protein_sequence': protein_sequence,
            'ph': ph,
            'd20_saturation': saturation,
            'temp': temperature,
            'file_type': file_type
        }
    include_exclude = None
    if PEPTIDE_LIST_PATH:
        include_exclude = args.include_exclude
        if include_exclude == "include":
            include_exclude = True
        else:
            include_exclude = False
    parser_map = {
        "DynamX": _parse_dynamX,
        "HDXworkbench": _parse_HDXWorkbench,
        "BioPharma": _parse_biopharma,
        "Byos": _parse_byos,
        "HDExaminer": _parse_HDExaminer,
        "Custom": _parse_custom,
    }


    if file_type in parser_map:
        hxms = _conver_PFhxms_to_hxms(INPUT_PATH,protein_sequence,saturation,ph,temperature,RAW_SPEC_PATH,protein_name)
        if hxms:
            hxms_good = write_hxms_file(hxms, OUTPUT_PATH, flags, None, None, PEPTIDE_LIST_PATH, include_exclude, None, None,True)
            if not hxms_good:
                print("Error while generating file. Please check your inputs!")


if __name__ == "__main__":
    main()