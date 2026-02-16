
from pathlib import Path
from datetime import datetime
from typing import Dict, Optional, Any


class HxmsData:
    def __init__(self):
        self.proteins = []
        self.state = {}
        self.peptides = {}
        self.match = None
        self.metadata = {
            'protein_sequence': "",
            'protein_name': "",
            'protein_state': "",
            'temperature': "",
            'ph': "",
            'saturation': "",
            'file_type': ""
        }



def write_hxms_file(hxms_data: HxmsData, output_path: str, flags: Optional[Dict[str, Any]] = None,protein_state_rename:str = None, protein_name_rename:str = None,peptide_list_parsed=None,include_exclude_bool=True,replicate_ID: int = None,bimodel_group: str = None,over_ride_save=None, save_match=False, save_fine_match=False) -> bool:
    from Helper_Functions import _sort_hxms_data_cc, _get_avg_zero_uptake
    hxms_data = _sort_hxms_data_cc(hxms_data)

    metadata = hxms_data.metadata

    if flags:
        if not isinstance(flags['protein_sequence'], list):
            metadata['protein_sequence'] = flags.get('protein_sequence', "")
        metadata['protein_sequence'] = flags.get('protein_sequence', "")
        metadata['ph'] = flags.get('ph', "")
        metadata['saturation'] = flags.get('d20_saturation', "")
        metadata['temperature'] = flags.get('temp', "")

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
            if num_sequences == len(hxms_data.proteins):
                try:
                    metadata['protein_sequence'] = flags['protein_sequence'][index]
                except IndexError:
                    metadata['protein_sequence'] = flags['protein_sequence'][-1]
            elif num_sequences == 1:
                metadata['protein_sequence'] = flags['protein_sequence'][0]

            else:
                try:
                    metadata['protein_sequence'] = flags['protein_sequence'][index]
                except:
                    metadata['protein_sequence'] = flags['protein_sequence'][-1]
        elif isinstance(flags.get('protein_sequence'), str):
            metadata['protein_sequence'] = flags['protein_sequence']
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
                match_lines = []
                if protein in hxms_data.peptides and state in hxms_data.peptides[protein]:
                    for peptide_name in hxms_data.peptides[protein][state]:
                        if peptide_list_parsed is not None:
                            # peptide_name = f"{row['start']}-{peptide}-{modification}-{row['end']}-{replicate}"
                            start = peptide_name.split("~")[0]
                            end = peptide_name.split("~")[3]
                            key = start + "~" + end
                            if include_exclude_bool: #Include
                                if key not in peptide_list_parsed:
                                    continue
                            else: #Exclude
                                if key in peptide_list_parsed:
                                    continue

                        PTM_name = peptide_name.split("~")[2]
                        if PTM_name == "N/A":
                            PTM_name = "nan"
                        if PTM_name not in PTM_dic:
                            PTM_dic[PTM_name] = f"0000{len(PTM_dic)}"[-4:]


                        if ("0" not in hxms_data.peptides[protein][state][peptide_name]) and ("0.0" not in hxms_data.peptides[protein][state][peptide_name]):
                            zero_uptake = _get_avg_zero_uptake(hxms_data, protein, state, peptide_name)

                        else:
                            count = 0
                            zero_uptake = 0
                            if "0" in hxms_data.peptides[protein][state][peptide_name]:
                                for x in hxms_data.peptides[protein][state][peptide_name]["0"]:
                                    count += 1
                                    zero_uptake += x['uptake']
                            if "0.0" in hxms_data.peptides[protein][state][peptide_name]:
                                for x in hxms_data.peptides[protein][state][peptide_name]["0.0"]:
                                    count += 1
                                    zero_uptake += x['uptake']
                            if count != 0:
                                zero_uptake = zero_uptake/count

                        for exposure_key in hxms_data.peptides[protein][state][peptide_name]:
                            for tp in hxms_data.peptides[protein][state][peptide_name][exposure_key]:
                                index += 1
                                uptake = tp['uptake']
                                if ((exposure_key != "0") and (exposure_key != "0.0")) and zero_uptake is not None:
                                    uptake = tp['uptake'] - zero_uptake
                                elif (exposure_key == "0") or (exposure_key == "0.0"):
                                    # Zero-time point uptake is always reported as 0
                                    uptake = 0
                                else:
                                    continue
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
                                    f"{float(tp['time']):<16.6e}"
                                    f"{uptake:<9.2f}"
                                    f"{tp['envelope']}\n"
                                )
                                tp_lines.append(line)
                                
                                if save_match:
                                    mact_content = tp['raw_ms_fine_structure'] if save_fine_match else tp['match']
                                    match_line = (
                                        f"{'MATCH':<12}"
                                        f"{index:<8}"
                                        f"{tp['score']:<8}"
                                        f"{tp['RT']:<8}"
                                        f"{tp['charge_state']:<5}"
                                        f"{tp['mono_mz']:<16}"
                                        #f"{tp['match']}\n"
                                        f"{mact_content}\n"
                                    )
                                    match_lines.append(match_line)

                col_title_PTM = 'TITLE_PTM   PTM_ID  CONTENT\n'
                col_title_match = 'TITLE_MATCH TP_ID   CONF    RT(min) Z    MONO_M          M/Z\n'
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
                        if save_match:
                            f.write(col_title_match)
                            f.writelines(match_lines)
                    if len(tp_lines) < 3:
                        return False


                else:
                    with open(output_path, 'w') as f:
                        f.write(header)
                        f.write(col_title_tp)
                        f.writelines(tp_lines)
                        f.write(col_title_PTM)
                        f.writelines(PTM_lines)
                        if save_match or save_fine_match:
                            f.write(col_title_match)
                            f.writelines(match_lines)
                    if len(tp_lines) < 3:
                        return False
    return True


def write_hxms_file_combined_test(hxms_data: HxmsData, output_path: str, ) -> bool:

    metadata = hxms_data.metadata
    files_written = 0
    for protein in hxms_data.proteins:
        if protein in hxms_data.state:
            for state in hxms_data.state[protein]:

                metadata['protein_name'] = protein
                metadata['protein_state'] = state
                PTM_dic = {"nan": "0000"}
                sortable_tp_lines = []
                sortable_match_lines = []
                has_match_data = False

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

                if protein in hxms_data.peptides and state in hxms_data.peptides[protein]:
                    index = 0
                    for peptide_name in hxms_data.peptides[protein][state]:

                        PTM_name = peptide_name.split("~")[2]
                        if PTM_name == "N/A":
                            PTM_name = "nan"
                        if PTM_name not in PTM_dic:
                            PTM_dic[PTM_name] = f"0000{len(PTM_dic)}"[-4:]

                        for exposure_key in hxms_data.peptides[protein][state][peptide_name]:

                            for tp in hxms_data.peptides[protein][state][peptide_name][exposure_key]:
                                index += 1
                                uptake = tp['uptake']

                                time_val = tp['time']
                                time_str = 'inf' if time_val == float('inf') else f"{time_val:.6e}"
                                line = (
                                    f"{'TP':<12}"
                                    f"{index:<8}"
                                    f"{tp['mod']:<6}"
                                    f"{tp['start']:<7}"
                                    f"{tp['end']:<7}"
                                    f"{tp['replicate']:<5}"
                                    f"{PTM_dic[PTM_name]:<8}"
                                    f"{float(tp['time']):<16.6e}"
                                    f"{uptake:<9.2f}"
                                    f"{tp['envelope']}\n"
                                )

                                sort_key = (
                                    tp['start'],
                                    tp['end'],
                                    tp['replicate'],
                                    time_val if time_val != 'inf' else float('inf')
                                )

                                sortable_tp_lines.append((sort_key, line))
                                
                                if 'match' in tp and tp['match']:
                                    has_match_data = True
                                    rt_val = tp.get('RT', '')
                                    charge_val = tp.get('charge_state', '')
                                    mono_mz_val = tp.get('mono_mz', '')
                                    score = tp.get('score', '')
                                    match_content = tp.get('raw_ms_fine_structure', tp['match'])
                                    
                                    match_line = (
                                        f"{'MATCH':<12}"
                                        f"{index:<8}"
                                        f"{score:<8}"
                                        f"{rt_val:<8}"
                                        f"{charge_val:<5}"
                                        f"{mono_mz_val:<16}"
                                        f"{match_content}\n"
                                    )
                                    sortable_match_lines.append((sort_key, match_line))

                sortable_tp_lines.sort(key=lambda item: item[0])
                final_tp_lines = [item[1] for item in sortable_tp_lines]
                
                if has_match_data:
                    sortable_match_lines.sort(key=lambda item: item[0])
                    final_match_lines = [item[1] for item in sortable_match_lines]

                col_title_PTM = 'TITLE_PTM  PTM_ID  CONTENT\n'
                col_title_match = 'TITLE_MATCH TP_ID   CONF    RT(min) Z    MONO_M          M/Z\n'
                PTM_lines = []
                for PTM in PTM_dic:
                    line = f"{'PTM':<12}{PTM_dic[PTM]:<8}{PTM}\n"
                    PTM_lines.append(line)
                current_output_path = output_path

                with open(current_output_path, 'w') as f:
                    f.write(header)
                    f.write(col_title_tp)
                    f.writelines(final_tp_lines)
                    f.write(col_title_PTM)
                    f.writelines(PTM_lines)
                    if has_match_data:
                        f.write(col_title_match)
                        f.writelines(final_match_lines)

                files_written += 1
                if len(final_tp_lines) < 3:
                    return False

    return files_written > 0