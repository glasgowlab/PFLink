
from pathlib import Path
from datetime import datetime
from typing import Dict, Optional, Any


class HxmsData:
    """
    A class to store and structure data from HDX-MS experiments.
    """
    def __init__(self):
        # List of proteins found in the data
        self.proteins = []
        # Dictionary to store states for each protein
        self.state = {}
        # Dictionary to store peptide data, nested by protein and state
        self.peptides = {}
        # Dictionary for metadata, to be populated from the file or flags
        self.metadata = {
            'protein_sequence': "",
            'protein_name': "",
            'protein_state': "",
            'temperature': "",
            'ph': "",
            'saturation': ""
        }



def write_hxms_file(hxms_data: HxmsData, output_path: str, flags: Optional[Dict[str, Any]] = None,protein_state_rename:str = None, protein_name_rename:str = None,peptide_list_parsed=None,include_exclude_bool=True,replicate_ID: int = None,bimodel_group: str = None,over_ride_save=None) -> bool:
    """
    Writes the data from the HxmsData object to one or more .hxms files,
    based on the protein and state.
    """
    from Helper_Functions import _sort_hxms_data, _get_avg_zero_uptake
    hxms_data = _sort_hxms_data(hxms_data)

    metadata = hxms_data.metadata

    # Use flags to pre-populate metadata if available
    if flags:
        if not isinstance(flags['protein_sequence'], list):
            metadata['protein_sequence'] = flags.get('protein_sequence', "")
        metadata['protein_sequence'] = flags.get('protein_sequence', "")
        metadata['ph'] = flags.get('ph', "")
        metadata['saturation'] = flags.get('d20_saturation', "")
        metadata['temperature'] = flags.get('temp', "")

    # Logic to handle multiple proteins/states from flags
    protein_states_to_process = []
    if flags and flags.get('protein_name_states'):
        protein_name_states = flags['protein_name_states']
        if isinstance(protein_name_states, list):
            protein_states_to_process = [tuple(pair) for pair in protein_name_states]
        elif isinstance(protein_name_states, (list, tuple)) and len(protein_name_states) == 2:
            protein_states_to_process = [tuple(protein_name_states)]
        elif isinstance(protein_name_states, str):
            parts = protein_name_states.split(" ")
            if len(parts) >= 2:
                protein_states_to_process = [(parts[0], parts[1])]
            else:
                print(f"Warning: protein_name_states string '{protein_name_states}' could not be parsed.")
                protein_states_to_process = []

    num_sequences = len(flags.get('protein_sequence', []))
    for index, protein in enumerate(hxms_data.proteins):


        if isinstance(flags.get('protein_sequence'), list):
            # Scenario 1: Number of sequences matches the number of proteins (1-to-1 mapping)
            if num_sequences == len(hxms_data.proteins):
                try:
                    # Use the current index to map a sequence to the current protein
                    metadata['protein_sequence'] = flags['protein_sequence'][index]
                except IndexError:
                    # Fallback to the last one (safety, though logic should prevent this)
                    metadata['protein_sequence'] = flags['protein_sequence'][-1]

            # Scenario 2: Only one sequence is provided (use it for all proteins)
            elif num_sequences == 1:
                # Use the first (and only) sequence provided
                metadata['protein_sequence'] = flags['protein_sequence'][0]

            # Scenario 3: Otherwise (e.g., 2 sequences for 3 proteins, or 0)
            else:
                try:
                    metadata['protein_sequence'] = flags['protein_sequence'][index]
                except:
                    metadata['protein_sequence'] = flags['protein_sequence'][-1]


        # If protein_sequence is a single string (not a list), use it for all proteins.
        elif isinstance(flags.get('protein_sequence'), str):
            metadata['protein_sequence'] = flags['protein_sequence']

        # If no protein_sequence flag was set or it's an unsupported type, use an empty string.
        else:
            metadata['protein_sequence'] = ""

        if protein in hxms_data.state:
            for state in hxms_data.state[protein]:
                if protein_states_to_process and (protein, state) not in protein_states_to_process:
                    continue

                metadata['protein_name'] = protein
                metadata['protein_state'] = state

                if protein_state_rename is not None and protein_state_rename != "":
                    metadata['protein_state'] = protein_state_rename
                if protein_name_rename is not None and protein_name_rename != "":
                    metadata['protein_name'] = protein_name_rename

                PTM_dic = {"nan": "0000"}

                header = (
                    #f"HEADER    HX/MS DATA FORMAT v1.0 {datetime.now().strftime('%Y-%m-%d')}\n"
                    #f"HEADER    Hydrogen Exchange Mass Spectrometry Data\n"
                    #f"REMARK100\n"
                    #f"REMARK100 DOI:\n"
                    f"HEADER       HX/MS DATA FORMAT v1.0 {datetime.now().strftime('%Y-%m-%d')}\n"
                    f"HEADER       Hydrogen Exchange Mass Spectrometry Data\n"
                    f"METADATA     PROTEIN_SEQUENCE      {metadata['protein_sequence']}\n"
                    f"METADATA     PROTEIN_NAME          {metadata['protein_name']}\n"
                    f"METADATA     PROTEIN_STATE         {metadata['protein_state']}\n"
                    f"METADATA     TEMPERATURE(K)        {metadata['temperature']}\n"
                    f"METADATA     pH(READ)              {metadata['ph']}\n"
                    f"METADATA     D2O_SATURATION        {metadata['saturation']}\n"
                )
                col_title_tp = (
                    f"{'TITLE_TP':<12}"
                    f"{'INDEX':<8}"
                    f"{'MOD':<6}"
                    f"{'START':<7}"
                    f"{'END':<7}"
                    f"{'REP':<5}"
                    f"{'PTM_ID':<8}"
                    f"{'TIME(Sec)':<16}"
                    f"{'UPTAKE':<9}"
                    f"ENVELOPE\n"
                )

                tp_lines = []
                index = 0
                if protein in hxms_data.peptides and state in hxms_data.peptides[protein]:
                    for peptide_name in hxms_data.peptides[protein][state]:
                        if peptide_list_parsed is not None:
                            # peptide_name = f"{row['start']}-{peptide}-{modification}-{row['end']}-{replicate}"
                            start = peptide_name.split("-")[0]
                            end = peptide_name.split("-")[3]
                            key = start + "-" + end
                            if include_exclude_bool: #Include
                                if key not in peptide_list_parsed:
                                    continue
                            else: #Exclude
                                if key in peptide_list_parsed:
                                    continue

                        PTM_name = peptide_name.split("-")[2]
                        if PTM_name == "N/A":
                            PTM_name = "nan"
                        if PTM_name not in PTM_dic:
                            PTM_dic[PTM_name] = f"0000{len(PTM_dic) - 1}"[-4:]

                        # Check for a timepoint "0" and skip if it doesn't exist.
                        if "0" not in hxms_data.peptides[protein][state][peptide_name]:
                            zero_uptake = _get_avg_zero_uptake(hxms_data, protein, state, peptide_name)

                        else:
                            zero_uptake = hxms_data.peptides[protein][state][peptide_name]["0"]['uptake']

                        for exposure_key in hxms_data.peptides[protein][state][peptide_name]:
                            index += 1
                            tp = hxms_data.peptides[protein][state][peptide_name][exposure_key]
                            uptake = tp['uptake']
                            if exposure_key != "0" and zero_uptake is not None:
                                uptake = tp['uptake'] - zero_uptake
                            elif exposure_key == "0":
                                # Zero-time point uptake is always reported as 0
                                uptake = 0
                            if replicate_ID:
                                tp['replicate'] = replicate_ID
                            if bimodel_group:
                                tp['mod'] = replicate_ID


                            line = (
                                f"{'TP':<12}"
                                f"{index:<8}"
                                f"{tp['mod']:<6}"
                                f"{tp['start']:<7}"
                                f"{tp['end']:<7}"
                                f"{tp['replicate']:<5}"
                                f"{PTM_dic[PTM_name]:<8}"
                                f"{tp['time']:<16.6e}"
                                f"{uptake:<9.2f}"
                                f"{tp['envelope']}\n"
                            )
                            tp_lines.append(line)

                col_title_PTM = 'TITLE_PTM  PTM_ID  CONTENT\n'
                PTM_lines = []
                for PTM in PTM_dic:
                    line = f"{'PTM':<12}{PTM_dic[PTM]:<8}{PTM}\n"
                    PTM_lines.append(line)
                if over_ride_save:
                    with open(f'{output_path}/{protein}_{state}.hxms', 'w') as f:
                        f.write(header)
                        f.write(col_title_tp)
                        f.writelines(tp_lines)
                        f.write(col_title_PTM)
                        f.writelines(PTM_lines)
                    if len(tp_lines) < 3:
                        return False


                else:
                    with open(output_path, 'w') as f:
                        f.write(header)
                        f.write(col_title_tp)
                        f.writelines(tp_lines)
                        f.write(col_title_PTM)
                        f.writelines(PTM_lines)
                    if len(tp_lines) < 3:
                        return False
    return True


def write_hxms_file_combined_test(hxms_data: HxmsData, output_path: str, ) -> bool:
    """
    Writes the data from the HxmsData object to one or more .hxms files,
    sorting the output lines by START -> END -> REPLICATE -> TIME.
    """
    # from Helper_Functions import _sort_hxms_data, _get_avg_zero_uptake # Assuming available

    metadata = hxms_data.metadata
    files_written = 0

    # --------------------------------------------------------------------------
    # Outer loops for Protein and State (defines the file content and name)
    # --------------------------------------------------------------------------
    for protein in hxms_data.proteins:
        if protein in hxms_data.state:
            for state in hxms_data.state[protein]:

                # --- File Setup ---
                metadata['protein_name'] = protein
                metadata['protein_state'] = state
                PTM_dic = {"nan": "0000"}

                # List to store (sort_key, formatted_line_string) tuples
                sortable_tp_lines = []

                header = (
                    f"HEADER       HX/MS DATA FORMAT v1.0 {datetime.now().strftime('%Y-%m-%d')}\n"
                    f"HEADER       Hydrogen Exchange Mass Spectrometry Data\n"
                    f"METADATA     PROTEIN_SEQUENCE      {metadata['protein_sequence']}\n"
                    f"METADATA     PROTEIN_NAME          {metadata['protein_name']}\n"
                    f"METADATA     PROTEIN_STATE         {metadata['protein_state']}\n"
                    f"METADATA     TEMPERATURE(K)        {metadata['temperature']}\n"
                    f"METADATA     pH(READ)              {metadata['ph']}\n"
                    f"METADATA     D2O_SATURATION        {metadata['saturation']}\n"
                )
                col_title_tp = (
                    f"{'TITLE_TP':<12}"
                    f"{'INDEX':<8}"
                    f"{'MOD':<6}"
                    f"{'START':<7}"
                    f"{'END':<7}"
                    f"{'REP':<5}"
                    f"{'PTM_ID':<8}"
                    f"{'TIME(Sec)':<16}"
                    f"{'UPTAKE':<9}"
                    f"ENVELOPE\n"
                )

                # --- Data Gathering and Formatting ---
                if protein in hxms_data.peptides and state in hxms_data.peptides[protein]:
                    # Initialize index outside, as it's a global counter for TP lines
                    index = 0

                    # Iterate through the data structure in its current order
                    for peptide_name in hxms_data.peptides[protein][state]:

                        # PTM Logic (remains the same)
                        PTM_name = peptide_name.split("-")[2]
                        if PTM_name == "N/A":
                            PTM_name = "nan"
                        if PTM_name not in PTM_dic:
                            PTM_dic[PTM_name] = f"{(len(PTM_dic) - 1):04d}"  # Cleaner format

                        for exposure_key in hxms_data.peptides[protein][state][peptide_name]:
                            for tp in hxms_data.peptides[protein][state][peptide_name][exposure_key]:
                                index += 1
                                uptake = tp['uptake']

                                # Handle 'inf' time formatting separately from scientific notation
                                time_val = tp['time']
                                time_str = 'inf' if time_val == float('inf') else f"{time_val:.6e}"

                                # Create the formatted line string
                                line = (
                                    f"{'TP':<12}"
                                    f"{index:<8}"
                                    f"{tp['mod']:<6}"
                                    f"{tp['start']:<7}"
                                    f"{tp['end']:<7}"
                                    f"{tp['replicate']:<5}"
                                    f"{PTM_dic[PTM_name]:<8}"
                                    f"{time_str:<16}"
                                    f"{uptake:<9.2f}"
                                    f"{tp['envelope']}\n"
                                )

                                # Create the sort key tuple: (START, END, REPLICATE, TIME)
                                sort_key = (
                                    tp['start'],
                                    tp['end'],
                                    tp['replicate'],
                                    # Use float('inf') for sorting if time is 'inf'
                                    time_val if time_val != 'inf' else float('inf')
                                )

                                # Store the key and the line together
                                sortable_tp_lines.append((sort_key, line))

                # --- Sorting and Final Line Assembly ---

                # Sort the list based on the sort_key (index 0 of the tuple)
                sortable_tp_lines.sort(key=lambda item: item[0])

                # Extract only the final formatted line strings
                final_tp_lines = [item[1] for item in sortable_tp_lines]

                # --- PTM Lines ---
                col_title_PTM = 'TITLE_PTM  PTM_ID  CONTENT\n'
                PTM_lines = []
                for PTM in PTM_dic:
                    line = f"{'PTM':<12}{PTM_dic[PTM]:<8}{PTM}\n"
                    PTM_lines.append(line)

                # --- File Writing ---
                # Consider generating a file name per protein/state if the intention is to write multiple files
                current_output_path = output_path

                with open(current_output_path, 'w') as f:
                    f.write(header)
                    f.write(col_title_tp)
                    f.writelines(final_tp_lines)  # Write the sorted lines
                    f.write(col_title_PTM)
                    f.writelines(PTM_lines)

                files_written += 1
                if len(final_tp_lines) < 3:
                    return False

    return files_written > 0