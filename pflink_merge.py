
import argparse
import sys
from rich_argparse import RichHelpFormatter
import os

from HXMS_IO import HxmsData
from Helper_Functions import combine_hxms_data
from Parsers import validate_and_parse_hxms_file
import tempfile
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
    parser.add_argument("--input_hxms_path1", type=str, required=True,
                        help="Path to the file1 containing the.hmxs data")
    parser.add_argument("--input_hxms_path2", type=str, required=True,
                        help="Path to the file2 containing the.hmxs data")
    parser.add_argument("--output_hxms_path", type=str, required=True,
                        help="Path to the output directory where the final .hxms files will be saved.")
    
    args = parser.parse_args()
    # Construct the full paths from the parsed arguments
    INPUT_PATH1 = args.input_hxms_path1
    INPUT_PATH2 = args.input_hxms_path1
    OUTPUT_PATH = args.output_hxms_path
    # Path Validations for input files
    if not os.path.isfile(INPUT_PATH1):
        print(f"{RED}Error: Your input path '{INPUT_PATH1}' does not exist or is not a file! Please correct it and run again.{RESET}")
        sys.exit(1)
    if not os.path.isfile(INPUT_PATH2):
        print(f"{RED}Error: Your input path '{INPUT_PATH2}' does not exist or is not a file! Please correct it and run again.{RESET}")
        sys.exit(1)
    if not os.path.isdir(OUTPUT_PATH):
        print(f"{RED}Error: Your output path '{OUTPUT_PATH}' does not exist or is not a directory! Please correct it and run again.{RESET}")
        sys.exit(1)

    hxms_data_1 = validate_and_parse_hxms_file(INPUT_PATH1)
    hxms_data_2 = validate_and_parse_hxms_file(INPUT_PATH2)
    if (not hxms_data_1[0]) or (not hxms_data_2[0]):
        print("One of your hxms files are wrong. Please fix them before trying to merge.")
        exit()

    hxms_1_hxms_object: HxmsData = hxms_data_1[2]
    hxms_2_hxms_object: HxmsData = hxms_data_2[2]
    
    if hxms_1_hxms_object.proteins != hxms_2_hxms_object.proteins:
        print("Your protein names do no match. Please sure you are using the same protein and try again.")
        exit()
    if hxms_1_hxms_object.state != hxms_2_hxms_object.state:
        print("Your protein states do no match. Please sure you are using the same state and try again.")
        exit()
    if hxms_1_hxms_object.metadata['protein_sequence'] != hxms_2_hxms_object.metadata['protein_sequence']:
        print("Your protein sequences do no match. Please sure you are using the same sequence and try again.")
        exit()
    if hxms_1_hxms_object.metadata['temperature'] != hxms_2_hxms_object.metadata['temperature']:
        print("Your temperatures do no match. Please sure you are using the same temperature and try again.")
        exit()
    if hxms_1_hxms_object.metadata['ph'] != hxms_2_hxms_object.metadata['ph']:
        print("Your phs do no match. Please sure you are using the same phs and try again.", duration=30)
        exit()
    if hxms_1_hxms_object.metadata['saturation'] != hxms_2_hxms_object.metadata['saturation']:
        print("Your saturations do no match. Please sure you are using the same saturations and try again.")
        exit()

    out_hxms = combine_hxms_data(hxms_1_hxms_object, hxms_data_1[-1], hxms_2_hxms_object, hxms_data_2[-1])

    from HXMS_IO import write_hxms_file_combined_test
    write_hxms_file_combined_test(out_hxms, f"{hxms_1_hxms_object.metadata['protein_name']}_{hxms_1_hxms_object.metadata['protein_state']}_merged.hxms")



if __name__ == "__main__":
    main()