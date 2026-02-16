
import os
import pandas as pd
import ast
import csv
from HXMS_IO import HxmsData
from Helper_Functions import validate_protein_sequence, _convert_dformula_to_dict, _create_peptide_chemical_formula, _subtract_chemical_formula_dicts

def _parse_peptide_list(file_path: str) -> list:
    if not os.path.exists(file_path) or not os.path.isfile(file_path):
        return None

    try:
        peptide_df = pd.read_csv(file_path)
        peptide_df.columns = peptide_df.columns.str.lower()
        if 'start' not in peptide_df.columns or 'end' not in peptide_df.columns:
            return None

        output_peptide_list = []
        for _, row in peptide_df.iterrows():
            start = str(row['start']).strip()
            end = str(row['end']).strip()
            output_peptide_list.append(f"{start}~{end}")

        return output_peptide_list
    except Exception as e:
        return None



class FlagsParser:
    def __init__(self, file_path: str):
        self.file_path = file_path
        self.flags = {}

    def _read_flags(self) -> bool:
        if not os.path.exists(self.file_path) or not os.path.isfile(self.file_path):
            return False
        try:
            with open(self.file_path, 'r') as flags_file:
                for line in flags_file:
                    if ":" in line:
                        key = line.split(":")[0]
                        value = ast.literal_eval(line.split(":")[1])
                        self.flags[key] = value
            return True
        except:
            return False

    def _validate_flags_file(self) -> [bool,str]:
        return_data = [True,""]
        if len(self.flags) == 0:
            return_data = [False, "You must have a protein state and name"]
        for key in self.flags:
            if key == 'protein_name_states':
                if not isinstance(self.flags[key], list):
                    return_data = [False, "protein_name_states must be a list"]
                    break
                if len(self.flags[key]) == 0:
                    return_data = [False, "You must have a protein state and name"]
                    break
                for item in self.flags[key]:
                    if len(item) != 2:
                        return_data = [False, "You must have a protein state and name for each entry"]
                        break
                    if (not isinstance(item[0], str)) or (not isinstance(item[1], str)):
                        return_data = [False, "The protein name and state must be a string"]
                        break

            elif key =='protein_sequence':
                if not isinstance(self.flags[key], list):
                    return_data = [False, "Protein sequences must be a list"]
                    break
                if len(self.flags[key]) == 0:
                    return_data = [False, "You must have at least one protein sequence"]
                    break
                for item in self.flags[key]:
                    if not isinstance(item, str):
                        return_data = [False, "The protein sequences must be a string"]
                        break
                    if not self._validate_protein_sequence(item):
                        return_data = [False, "The protein sequences must contain only the 20 native AAs"]
                        break
            elif key == 'ph':
                if not isinstance(self.flags[key], float) and not isinstance(self.flags[key], int):
                    return_data = [False, "The pH must be a float"]
                    break
                if self.flags[key] > 14 or self.flags[key] < 0:
                    return_data = [False, "The pH must be between 0 and 14"]
                    break

            elif key == 'd20_saturation':
                if not isinstance(self.flags[key], float) and not isinstance(self.flags[key], int):
                    return_data = [False, "D20_saturation must be a float"]
                    break
                if self.flags[key] > 1.0 or self.flags[key] < 0.0:
                    return_data = [False, "D20_saturation must be between 0 and 1"]
                    break

            elif key == 'temp':
                if not isinstance(self.flags[key], float) and not isinstance(self.flags[key], int):
                    return_data = [False, "Temp must be a float or int"]
                    break

            elif key == 'file_type':
                if not isinstance(self.flags[key], str):
                    return_data = [False, "file_type must be string"]
                    break
                if self.flags[key] not in ["DynamX", "HDXworkbench", "HDExaminer", "BioPharma", "Custom"]:
                    return_data = [False, "file_type must be a supported filed type"]
                break
        if not return_data[0]:
            raise SyntaxError(f"Invalid format {return_data[1]}")

        return return_data
    def _validate_protein_sequence(self, sequence: str) -> bool:
        allowed_chars = "ACDEFGHIKLMNPQRSTVWY"
        if not all(c.upper() in allowed_chars for c in sequence):
            return False
        return True

    def parse(self) -> dict:
        if not self._read_flags():
            return None
        self._validate_flags_file()
        return self.flags

def _parse_HDXWorkbench(file_path: str, flags_info: None) -> HxmsData:
    """
    Parses an HDXworkbench CSV file into a HxmsData object.
    """
    hx_data = HxmsData()
    if not os.path.exists(file_path) or not os.path.isfile(file_path):
        return None

    hdxworkbench_important_headers = [
        "peptide", "charge", "start", "end", "numExHydrogens", "experiment",
        "sample", "timepoint", "replicate", "discarded_replicate",
        "percentd_replicate", "centroid",'dformula'
    ]
    metadata = {}
    header_row_index = None
    with open(file_path, 'r', newline='', encoding='utf-8') as file:
        reader = csv.reader(file)
        for i, row in enumerate(reader):
            if row and len(row) > 1 and row[0].strip().lower() == 'peptide' and row[1].strip().lower() == 'charge':
                header_row_index = i
                break
            if row and len(row) >= 2 and row[0].strip():
                key = row[0].strip()
                value = row[1].strip()
                metadata[key] = value

    if flags_info is None:
        hx_data.metadata['file_type'] = "HDXworkbench"
        for key in metadata:
            if key.lower() == "ph":
                hx_data.metadata["ph"] = str(metadata[key])
            if key.lower() == "deuterium solution concentration":
                hx_data.metadata["saturation"] = str(metadata[key])
            if key.lower() == "temperature":
                hx_data.metadata["temperature"] = str(float(metadata[key]) + 273.15)
            if key.lower() == "experiment protein sequence":
                hx_data.metadata["protein_sequence"] = str(metadata[key])

    if header_row_index is not None:
        try:
            hdxworkbench_df = pd.read_csv(file_path, skiprows=header_row_index)
        except Exception as e:
            print(f"Error reading file: {e}")
            return None

        for header in hdxworkbench_important_headers:
            if header not in hdxworkbench_df.columns:
                return None

        for _, row in hdxworkbench_df.iterrows():
            state = str(row['sample'])
            if len(state.split(" ")) > 1:
                protein = state.split(" ")[0]
                state = state.split(" ")[1]
            elif state == "nan":
                continue

            modification = ""
            peptide = str(row['peptide'])
            if "[" in peptide:
                peptide = row['peptide'].split("[")[0]
                modification = row['peptide'].split("[")[1].split("]")[0]

            dformula_created = _create_peptide_chemical_formula(peptide)
            dformula = _convert_dformula_to_dict(row['dformula'])
            modification = _subtract_chemical_formula_dicts(dformula_created,dformula)
            if modification == "":
                modification = "N/A"
            modification = modification.replace("-", "_")
            envelope = ""
            replicate = int(row["replicate"]) - 1
            if "isotope" in hdxworkbench_df.columns:
                envelope = row["isotope"]
                envelope = envelope.replace("NA","0.00")
                envelope = envelope.replace(":",",")
            peptide_name = f"{row['start']}~{peptide}~{modification}~{row['end']}~{replicate}"
            protein = "protein"
            if protein not in hx_data.proteins:
                hx_data.proteins.append(protein)
                hx_data.state[protein] = []
                hx_data.peptides[protein] = {}
            if state not in hx_data.state[protein]:
                hx_data.state[protein].append(state)
                hx_data.peptides[protein][state] = {}
            if peptide_name not in hx_data.peptides[protein][state]:
                hx_data.peptides[protein][state][peptide_name] = {}

            if str(row["timepoint"]) != "NA" and str(row["timepoint"]) != "nan":
                exposure = str(row["timepoint"])[0:-1]
                if exposure == "999999":
                    exposure = "inf"
                if str(row["discarded_replicate"]) == "True":
                    continue

                if exposure not in hx_data.peptides[protein][state][peptide_name]:

                    hx_data.peptides[protein][state][peptide_name][exposure] = [{
                        "uptake": float(row['centroid']) * float(row['charge']),
                        "start": int(row['start']),
                        "end": int(row['end']),
                        "sequence": str(row['peptide']),
                        "PTM": modification,
                        "time": float(exposure),
                        "mod": "A",
                        "envelope": envelope,
                        "replicate": int(row["replicate"]) - 1
                    }]
                else:
                    hx_data.peptides[protein][state][peptide_name][exposure].append({
                        "uptake": float(row['centroid']) * float(row['charge']),
                        "start": int(row['start']),
                        "end": int(row['end']),
                        "sequence": str(row['peptide']),
                        "PTM": modification,
                        "time": float(exposure),
                        "mod": "A",
                        "envelope": envelope,
                        "replicate": int(row["replicate"]) - 1
                    })
    return hx_data

def _sort_hxms_data_cc(hx_data: HxmsData) -> HxmsData:

    hx_data.proteins.sort()
    for protein in hx_data.state:
        hx_data.state[protein].sort()

    for protein in hx_data.peptides:
        for state in hx_data.peptides[protein]:

            # 'START-SEQ-PTM_ID-END'
            for _, time_data in hx_data.peptides[protein][state].items():
                time_data = {sorted(time_data.items(), key=lambda item:item[0])}

                for tp, tp_data in time_data.items():
                    tp_data = sorted(tp_data, key=lambda i: (i['replicate'],i['time']))


            hx_data.peptides[protein][state] = {
            sorted(
                hx_data.peptides[protein][state].items(),
                key=lambda item: (int(item[0].split('~')[0]),int(item[0].split('~')[-1]),item[0].split('~')[2]
                ))
            }


    return hx_data



def _parse_HDExaminer(file_path: str, flags_info: None) -> HxmsData:
    """
    Parses an HDExaminer CSV file into a HxmsData object.
    """
    hx_data = HxmsData()
    if not os.path.exists(file_path) or not os.path.isfile(file_path):
        return None

    header_row_index = None
    time_points = []
    with open(file_path, 'r', newline='', encoding='utf-8') as file:
        for index, line in enumerate(file):
            if index == 0:
                parts = line.replace("\n", "").split(",")
                for item in parts:
                    if item == "":
                        continue
                    else:
                        if item == "Full-D":
                            time_points.append("inf")
                        else:
                            try:
                                float(item[0:-1])
                                time_points.append(item[0:-1])
                            except:
                                pass
            if index > 0:
                header_row_index = 1
                break
    if header_row_index is not None:
        try:
            hdx_df = pd.read_csv(file_path, skiprows=1)
        except Exception as e:
            print(f"Error reading file: {e}")
            return None

        hdx_df.columns = hdx_df.columns.str.strip()

        required_headers = ["State", "Protein", "Sequence", "Start", "End", "Charge"]
        for header in required_headers:
            if header not in hdx_df.columns:
                print(f"Missing required header: {header}")
                return None

        d_cols = [col for col in hdx_df.columns if (('#D' in col) and ('right' not in col))]

        for _, row in hdx_df.iterrows():
            used_rep = {}
            for item in time_points:
                used_rep[item]=0

            state = str(row['State']).strip()
            protein = str(row['Protein']).strip()
            if not protein or protein == "nan":
                protein = "protein"

            peptide = str(row['Sequence']).strip()
            start = int(row['Start'])
            end = int(row['End'])
            modification = "N/A"
            modification = modification.replace("-", "_")

            if protein not in hx_data.proteins:
                hx_data.proteins.append(protein)
                hx_data.state[protein] = []
                hx_data.peptides[protein] = {}
            if state not in hx_data.state[protein]:
                hx_data.state[protein].append(state)
                hx_data.peptides[protein][state] = {}


            for i, d_col in enumerate(d_cols):
                if pd.notna(row[d_col]):
                    time_point_str = time_points[i]
                    time_point = float(time_point_str) if time_point_str != "inf" else float('inf')
                    uptake = float(row[d_col])
                    used_rep[time_point_str]+=1
                    peptide_name = f"{row['Start']}~{row['Sequence']}~{modification}~{row['End']}~0"
                    #peptide_name = f"{row['Start']}-{row['Sequence']}-{modification}-{row['End']}-"
                    #peptide_name = peptide_name+str(used_rep[time_point_str]-1)
                    if peptide_name not in hx_data.peptides[protein][state]:
                        hx_data.peptides[protein][state][peptide_name] = {}

                    if "0.0" not in hx_data.peptides[protein][state][peptide_name]:
                        hx_data.peptides[protein][state][peptide_name]["0.0"] = [{

                            "uptake": 0,
                            "start": start,
                            "end": end,
                            "sequence": peptide,
                            "time": 0.0,
                            "envelope": "",
                            "PTM": modification,
                            "mod": "A",
                            #"replicate": used_rep[time_point_str]-1,
                            "replicate": 0
                        }]

                    if time_point not in hx_data.peptides[protein][state][peptide_name]:
                        hx_data.peptides[protein][state][peptide_name][str(time_point)] = [{
                            "uptake": uptake,
                            "start": start,
                            "end": end,
                            "sequence": peptide,
                            "time": time_point,
                            "envelope": "",
                            "PTM": modification,
                            "mod": "A",
                            #"replicate": used_rep[time_point_str]-1,
                            "replicate": 0
                        }]
                    else:
                        hx_data.peptides[protein][state][peptide_name][str(time_point)].append({
                            "uptake": uptake,
                            "start": start,
                            "end": end,
                            "sequence": peptide,
                            "time": time_point,
                            "PTM": modification,
                            "envelope": "",
                            "mod": "A",
                            #"replicate": used_rep[time_point_str]-1,
                            "replicate": 0
                        })

    return hx_data


def _parse_custom(file_path: str, flags_info: None) -> HxmsData:
    """
    Parses a custom CSV file into a HxmsData object.
    """
    hx_data = HxmsData()
    if not os.path.exists(file_path) or not os.path.isfile(file_path):
        return None
    custom_headers = ["Protein","State","Start","Stop","Replicate","Modification","Bimodal_group","Exposure","Uptake","Envelope","RT","Conf","Z","Mono_M","M/Z"]
    try:
        custom_df = pd.read_csv(file_path)
    except Exception as e:
        print(f"Error reading file: {e}")
        return None
    for header in custom_headers:
        if header not in custom_df.columns:
            return None
    for _, row in custom_df.iterrows():
        protein = str(row['Protein'])
        state = str(row['State'])
        if protein == "" or protein == "nan" or state == "" or state == "nan" or protein is None or state is None:
            continue
        modification = str(row["Modification"]) if pd.notna(row["Modification"]) else "N/A"
        modification = modification.replace("-", "_")
        replicate = str(row["Replicate"]) if (pd.notna(row["Replicate"]) and row["Replicate"] != "") else "0"
        replicate = int(float(replicate))
        peptide_name = f"{int(row['Start'])}~N/A~{modification}~{int(row['Stop'])}~{replicate}"
        rt = float(row["RT"]) if pd.notna(row["RT"]) else None
        z = float(row["Z"]) if pd.notna(row["Z"]) else None
        mono_mz = float(row["Mono_M"]) if pd.notna(row["Mono_M"]) else None
        score = float(row["Conf"]) if pd.notna(row["Conf"]) else None
        m_z = ""
        envelope = ""
        if row['Envelope'] != "":
            try:
                if ";" in row['Envelope']:
                    parts = row['Envelope'].replace("\n","").split(";")
                    if len(parts) > 1:
                        envelope = ",".join(parts)
                else:
                    envelope = str(row['Envelope'])
            except:
                pass
        if row['M/Z'] != "":
            try:
                if "|" in row['M/Z']:
                    parts = row['M/Z'].replace("\n","").split("|")
                    hx_data.match = True
                    if len(parts) > 1:
                        m_z = ",".join(parts)
                else:
                    m_z = str(row['M/Z'])
            except:
                pass
        if protein not in hx_data.proteins:
            hx_data.proteins.append(protein)
            hx_data.state[protein] = []
            hx_data.peptides[protein] = {}
        if state not in hx_data.state[protein]:
            hx_data.state[protein].append(state)
            hx_data.peptides[protein][state] = {}
        if peptide_name not in hx_data.peptides[protein][state]:
            hx_data.peptides[protein][state][peptide_name] = {}
        if str(row["Exposure"]) not in hx_data.peptides[protein][state][peptide_name]:

            hx_data.peptides[protein][state][peptide_name][str(row["Exposure"])] = [{
                "uptake": float(row["Uptake"]),
                "start": int(row['Start']),
                "end": int(row['Stop']),
                "PTM": modification,
                "time": float(row["Exposure"]),
                "mod": str(row["Bimodal_group"]) if str(row["Bimodal_group"]) != "" else "A",
                "envelope": envelope,
                "replicate": replicate,
                "match": m_z,
                "RT": rt,
                "score": score,
                "charge_state": z,
                "mono_mz": mono_mz,
                "raw_ms_fine_structure":m_z
                
            }]
        else:
            hx_data.peptides[protein][state][peptide_name][str(row["Exposure"])].append({
                "uptake": float(row["Uptake"]),
                "start": int(row['Start']),
                "end": int(row['Stop']),
                "PTM": modification,
                "time": float(row["Exposure"]),
                "mod": str(row["Bimodal_group"]) if str(row["Bimodal_group"]) != "" else "A",
                "envelope": envelope,
                "replicate": replicate,
                "match": m_z,
                "RT": rt,
                "score": score,
                "charge_state": z,
                "mono_mz": mono_mz,
                "raw_ms_fine_structure":m_z
            })

    return hx_data


def _parse_dynamX(file_path: str, flags_info: None) -> HxmsData:
    """
    Parses a DynamX CSV file into a HxmsData object.
    """
    hx_data = HxmsData()
    if not os.path.exists(file_path) or not os.path.isfile(file_path):
        return None
    dynamX_headers = ["Protein", "Start", "End", "Sequence", "Modification", "Fragment", "MaxUptake", "MHP", "State",
                      "Exposure", "Center", "Center SD", "Uptake", "Uptake SD", "RT", "RT SD"]
    try:
        dynamX_df = pd.read_csv(file_path)
    except Exception as e:
        print(f"Error reading file: {e}")
        return None
    for header in dynamX_headers:
        if header not in dynamX_df.columns:
            return None
    for _, row in dynamX_df.iterrows():
        protein = str(row['Protein'])
        state = str(row['State'])
        modification = str(row["Modification"]) if pd.notna(row["Modification"]) else "N/A"
        modification = modification.replace("-", "_")
        peptide_name = f"{row['Start']}~{row['Sequence']}~{modification}~{row['End']}~0"
        if protein not in hx_data.proteins:
            hx_data.proteins.append(protein)
            hx_data.state[protein] = []
            hx_data.peptides[protein] = {}
        if state not in hx_data.state[protein]:
            hx_data.state[protein].append(state)
            hx_data.peptides[protein][state] = {}
        if peptide_name not in hx_data.peptides[protein][state]:
            hx_data.peptides[protein][state][peptide_name] = {}

        if str(row["Exposure"]) not in  hx_data.peptides[protein][state][peptide_name]:
            hx_data.peptides[protein][state][peptide_name][str(row["Exposure"])] = [{
                "uptake": float(row["Uptake"]),
                "start": int(row['Start']),
                "end": int(row['End']),
                "sequence": str(row['Sequence']),
                "PTM": modification,
                "time": float(row["Exposure"]),
                "mod": "A",
                "envelope": "",
                "replicate": 0
            }]
        else:
            hx_data.peptides[protein][state][peptide_name][str(row["Exposure"])].append({
                "uptake": float(row["Uptake"]),
                "start": int(row['Start']),
                "end": int(row['End']),
                "sequence": str(row['Sequence']),
                "PTM": modification,
                "time": float(row["Exposure"]),
                "mod": "A",
                "envelope": "",
                "replicate": 0
            })
    return hx_data

def _parse_biopharma(file_path: str, flags_info: None) -> HxmsData:
    """
    Parses a BioPharma CSV file into a HxmsData object.
    """
    hx_data = HxmsData()
    if not os.path.exists(file_path) or not os.path.isfile(file_path):
        return None
    biopharam_headers = ["Retention Time", "Peptide", "Charge", "Residue From", "Residue To", "Ref", "0% Control",
                         "100% Control", "stdev"]
    try:
        biopharma_df = pd.read_csv(file_path, skiprows=1)
    except Exception as e:
        print(f"Error reading file: {e}")
        return None
    for header in biopharam_headers:
        if header not in biopharma_df.columns:
            return None
    time_points = []
    states = {}
    tp = {}
    all_indexes = []
    with open(file_path, 'r', newline='', encoding='utf-8') as file:
        for index, line in enumerate(file):
            if index == 0:
                line=line.replace("\n","").replace("\r","")
                parts = line.replace("\n", "").split(",")[1:]
                for sub_index,item in enumerate(parts):
                    if item == "":
                        continue
                    else:
                        if sub_index+1 not in states:

                            states[sub_index+1] = item
                        all_indexes.append(sub_index+1)
            elif index == 1:
                line=line.replace("\n","").replace("\r","")
                parts = line.replace("\n", "").split(",")
                for sub_index,item in enumerate(parts):

                    if item == "":
                        continue
                    failed = []
                    if sub_index in all_indexes:
                        try:
                            if states[sub_index] not in tp:
                                tp[states[sub_index]] = {}
                            tp[states[sub_index]][sub_index]=(float(item))
                        except:
                            failed.append(sub_index)
                            pass

                    for x in failed:
                        all_indexes.remove(x)
            else:
                break
    try:
        hdx_df = pd.read_csv(file_path, skiprows=1)
    except Exception as e:
        print(f"Error reading file: {e}")
        return None

    hdx_df.columns = hdx_df.columns.str.strip()
    protein = "protein"
    if protein not in hx_data.proteins:
        hx_data.proteins.append(protein)
        hx_data.state[protein] = []
        hx_data.peptides[protein] = {}

    for key in states:
        state = states[key]
        if state not in hx_data.state[protein]:
            hx_data.state[protein].append(state)
            hx_data.peptides[protein][state] = {}
    for _, row in hdx_df.iterrows():
        peptide = str(row['Peptide']).strip().split(" = ")[0].replace("-","_")
        try:
            start = int(row['Residue From'])
        except:
            continue
        end = int(row['Residue To'])
        modification = "N/A"
        modification = modification.replace("-", "_")
        peptide_name = f"{start}~{peptide}~{modification}~{end}~0"


        if pd.notna(row["100% Control"]):
            if pd.notna(row["0% Control"]):
                full_d = (float(row["100% Control"]) - float(row["0% Control"])) * float(row["Charge"])
                for key in tp:
                    if peptide_name not in hx_data.peptides[protein][key]:
                        hx_data.peptides[protein][key][peptide_name] = {}
                    hx_data.peptides[protein][key][peptide_name]["inf"] = {
                        "uptake": full_d,
                        "start": start,
                        "end": end,
                        "sequence": peptide,
                        "time": float('inf'),
                        "envelope": "",
                        "PTM": modification,
                        "mod": "A",
                        # "replicate": used_rep[time_point_str]-1,
                        "replicate": 0
                    }

        for i, d_col in enumerate(all_indexes):
            if pd.notna(row.iloc[d_col]):

                time_point = str(tp[states[d_col]][d_col])
                #time_point = float(time_point_str) if time_point_str != "inf" else float('inf')
                uptake = float(row.iloc[d_col])

                if peptide_name not in hx_data.peptides[protein][states[d_col]]:
                    hx_data.peptides[protein][states[d_col]][peptide_name] = {}
                if "0.0" not in hx_data.peptides[protein][states[d_col]][peptide_name]:
                    hx_data.peptides[protein][states[d_col]][peptide_name]["0.0"] = [{
                        "uptake": 0,
                        "start": start,
                        "end": end,
                        "sequence": peptide,
                        "time": 0.0,
                        "envelope": "",
                        "PTM": modification,
                        "mod": "A",
                        # "replicate": used_rep[time_point_str]-1,
                        "replicate": 0
                    }]

                if time_point not in hx_data.peptides[protein][states[d_col]][peptide_name]:
                    hx_data.peptides[protein][states[d_col]][peptide_name][time_point] = [{
                        "uptake": uptake,
                        "start": start,
                        "end": end,
                        "sequence": peptide,
                        "time": time_point,
                        "envelope": "",
                        "PTM": modification,
                        "mod": "A",
                        # "replicate": used_rep[time_point_str]-1,
                        "replicate": 0
                    }]
                else:

                    hx_data.peptides[protein][states[d_col]][peptide_name][time_point].append({
                        "uptake": uptake,
                        "start": start,
                        "end": end,
                        "sequence": peptide,
                        "time": time_point,
                        "PTM": modification,
                        "envelope": "",
                        "mod": "A",
                        # "replicate": used_rep[time_point_str]-1,
                        "replicate": 0
                    })
    return hx_data




def validate_and_parse_hxms_file(hxms_file_path:str)->[bool,str,HxmsData, dict]:
    if not os.path.exists(hxms_file_path):
        return [False,"File does not exist",None, None]
    if not os.path.isfile(hxms_file_path):
        return [False,"File is not a file",None, None]
    header_info = {}
    tp_info = {}
    ptm_info = {}
    match_info = {}
    skip_lines_TP = 0
    skip_lines_TP_done = False
    skip_lines_PTM = 0
    skip_lines_PTM_done = False
    skip_lines_MATCH = 0
    skip_lines_MATCH_done = False
    with open(hxms_file_path,"r") as hxms_file:
        for line in hxms_file:
            if "TITLE_TP" in line:
                if "TITLE_TP    INDEX   MOD   START  END    REP  PTM_ID  TIME(Sec)       UPTAKE" not in line:
                    return [False,"invalid TP header",None, None]
                skip_lines_TP_done = True
            elif not skip_lines_TP_done:
                skip_lines_TP +=1

                if "HEADER" in line:
                    continue
                if "METADATA" in line:
                    data = [data.replace("\n","") for data in line.split(" ") if data.replace("\n","") !=""]
                    try:
                        header_info[data[1]] = data[2]
                    except:
                        return [False,f"Incorrect header info for {data[1]} in {line}",None, None]
            elif skip_lines_TP_done:
                data = [data.replace("\n", "") for data in line.split(" ") if data.replace("\n", "") != ""]
                if "TP" == data[0]:
                    try:
                        envelope = ""
                        if len(data) == 10:
                            envelope = data[9]
                        if data[2]+"~"+data[3]+"~"+data[4]+"~"+data[5]+"~"+data[7] not in tp_info:
                            tp_info[data[2]+"~"+data[3]+"~"+data[4]+"~"+data[5]+"~"+data[7]] = [{

                                "INDEX":int(data[1]),
                                "MOD": data[2],
                                "START": int(data[3]),
                                "END": int(data[4]),
                                "REP": int(data[5]),
                                "PTM_ID": int(data[6]),
                                "TIME(Sec)": float(data[7]),
                                "UPTAKE": float(data[8]),
                                "ENVELOPE": envelope,

                            }]
                        else:
                            tp_info[data[2] + "~" + data[3] + "~" + data[4] + "~" + data[5] + "~" + data[7]].append({

                                "INDEX":int(data[1]),
                                "MOD": data[2],
                                "START": int(data[3]),
                                "END": int(data[4]),
                                "REP": int(data[5]),
                                "PTM_ID": int(data[6]),
                                "TIME(Sec)": float(data[7]),
                                "UPTAKE": float(data[8]),
                                "ENVELOPE": envelope,

                            })
                    except:
                        return [False, f"Incorrect TP info for {line}", None]

            if "TITLE_PTM" in line:
                if "TITLE_PTM   PTM_ID  CONTENT" not in line:
                    return [False,"invalid PTM header",None, None]
                skip_lines_PTM_done = True
            elif not skip_lines_PTM_done:
                skip_lines_PTM+=1
            elif skip_lines_PTM_done:
                data = [data.replace("\n", "") for data in line.split(" ") if data.replace("\n", "") != ""]
                if "PTM" == data[0]:
                    try:
                        ptm_info[int(data[1])] = data[2]
                    except:
                        return [False, f"Incorrect PTM info for {line}", None, None]
            
            if "TITLE_MATCH" in line:
                if "TITLE_MATCH TP_ID" not in line:
                    return [False,"invalid MATCH header",None, None]
                skip_lines_MATCH_done = True
            elif not skip_lines_MATCH_done:
                skip_lines_MATCH+=1
            elif skip_lines_MATCH_done:
                data = [data.replace("\n", "") for data in line.split(" ") if data.replace("\n", "") != ""]
                if "MATCH" == data[0]:
                    try:
                        match_content = line.split(None, 5)
                        print(match_content)
                        if len(match_content) >= 5:
                            tp_id = int(match_content[1])
                            score_val = match_content[2]
                            rt_val = match_content[3]
                            charge_val = match_content[4]
                            mono_mz_val = match_content[5].split(None, 1)[0] if len(match_content[5].split(None, 1)) > 0 else ""
                            match_data = match_content[5].split(None, 1)[1].strip() if len(match_content[5].split(None, 1)) > 1 else ""
                            
                            match_info[tp_id] = {
                                "RT": rt_val,
                                "score": score_val,
                                "charge_state": charge_val,
                                "mono_mz": mono_mz_val,
                                "match": match_data
                            }
                    except:
                        pass

    if "PROTEIN_SEQUENCE" not in header_info:
        return [False,"PROTEIN_SEQUENCE not in header",None, None]
    elif header_info["PROTEIN_SEQUENCE"] == "":
        return [False,"PROTEIN_SEQUENCE is empty",None, None]
    elif not validate_protein_sequence(header_info["PROTEIN_SEQUENCE"]):
        return [False, "PROTEIN_SEQUENCE is invalid", None, None]

    if "PROTEIN_NAME" not in header_info:
        return [False,"PROTEIN_NAME not in header",None, None]
    elif header_info["PROTEIN_NAME"] == "":
        return [False,"PROTEIN_NAME is empty",None, None]

    if "PROTEIN_STATE" not in header_info:
        return [False,"PROTEIN_STATE not in header",None, None]
    elif header_info["PROTEIN_STATE"] == "":
        return [False,"PROTEIN_STATE is empty",None, None]

    if "TEMPERATURE(K)" not in header_info:
        return [False,"TEMPERATURE(K) not in header",None, None]
    elif header_info["TEMPERATURE(K)"] == "":
        return [False,"TEMPERATURE(K) is empty",None, None]
    else:
        try:
            temp = float(header_info["TEMPERATURE(K)"])
            if temp > 350 or temp < 200:
                return [False, "TEMPERATURE(K) is not within normal range", None, None]
        except:
            return [False, "TEMPERATURE(K) cannot be converted to a float", None, None]


    if "pH(READ)" not in header_info:
        return [False,"pH(READ) not in header",None, None]
    elif header_info["pH(READ)"] == "":
        return [False,"pH(READ) is empty",None, None]
    else:
        try:
            ph = float(header_info["pH(READ)"])
            if ph > 14 or ph < 0:
                return [False, "pH(READ) is not within normal range", None, None]
        except:
            return [False, "pH(READ) cannot be converted to a float", None, None]

    if "D2O_SATURATION" not in header_info:
        return [False,"D2O_SATURATION not in header",None, None]
    elif header_info["D2O_SATURATION"] == "":
        return [False,"D2O_SATURATION is empty",None, None]
    else:
        try:
            sat = float(header_info["D2O_SATURATION"])
            if sat > 1.0 or sat < 0:
                return [False, "D2O_SATURATION is not within normal range", None, None]
        except:
            return [False, "D2O_SATURATION cannot be converted to a float", None, None]

    if len(ptm_info) == 0:
        return [False, "PTM dic is empty", None, None]
    if len(tp_info) == 0:
        return [False, "TP dic is empty", None, None]
    hx_data = HxmsData()
    hx_data.metadata['protein_sequence'] = header_info["PROTEIN_SEQUENCE"]
    hx_data.metadata['protein_name'] = header_info["PROTEIN_NAME"]
    hx_data.metadata['protein_state'] = header_info["PROTEIN_STATE"]
    hx_data.metadata['temperature'] = header_info["TEMPERATURE(K)"]
    hx_data.metadata['ph'] = header_info["pH(READ)"]
    hx_data.metadata['saturation'] = header_info["D2O_SATURATION"]

    protein = hx_data.metadata['protein_name']
    state = hx_data.metadata['protein_state']

    if protein not in hx_data.proteins:
        hx_data.proteins.append(protein)
        hx_data.state[protein] = []
        hx_data.peptides[protein] = {}
    if state not in hx_data.state[protein]:
        hx_data.state[protein].append(state)
        hx_data.peptides[protein][state] = {}

    for tp_key in tp_info:
        for tp in tp_info[tp_key]:
            try:
                ptm_name = ptm_info.get(tp['PTM_ID'], 'nan')
                peptide_name = f"{tp['START']}~N/A~{ptm_name}~{tp['END']}"
                if peptide_name not in hx_data.peptides[protein][state]:
                    hx_data.peptides[protein][state][peptide_name] = {}
                if str(tp['TIME(Sec)']) not in hx_data.peptides[protein][state][peptide_name]:
                    hx_data.peptides[protein][state][peptide_name][str(tp['TIME(Sec)'])] = []
                
                tp_data = {
                    "uptake": tp['UPTAKE'],
                    "start": tp['START'],
                    "end": tp['END'],
                    "PTM": ptm_name,
                    "time": tp["TIME(Sec)"],
                    "mod": tp["MOD"],
                    "envelope": tp["ENVELOPE"],
                    "replicate": tp["REP"]
                }
                
                if tp['INDEX'] in match_info:
                    tp_data["RT"] = match_info[tp['INDEX']]["RT"]
                    tp_data["charge_state"] = match_info[tp['INDEX']]["charge_state"]
                    tp_data["mono_mz"] = match_info[tp['INDEX']]["mono_mz"]
                    tp_data["match"] = match_info[tp['INDEX']]["match"]
                    tp_data["score"] = match_info[tp['INDEX']]["score"]
                
                hx_data.peptides[protein][state][peptide_name][str(tp['TIME(Sec)'])].append(tp_data)
            except:
                return [False, "TP could not be read or converted properly", None]
    return [True,"",hx_data,ptm_info]
