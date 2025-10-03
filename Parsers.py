
import os
import ast
import configparser
import pandas as pd
from typing import Dict, List, Optional, Union, Any
import csv
from HXMS_IO import HxmsData
from Helper_Functions import validate_protein_sequence, _convert_dformula_to_dict, _create_peptide_chemical_formula, _subtract_chemical_formula_dicts

def _parse_peptide_list(file_path: str) -> list:
    """
    Parses a CSV file and extracts a list of peptide ranges.

    Args:
        file_path (str): The path to the CSV file.

    Returns:
        list: A list of peptide ranges in "Start-End" format, or None if the file is invalid.
    """
    if not os.path.exists(file_path) or not os.path.isfile(file_path):
        return None

    try:
        peptide_df = pd.read_csv(file_path)
        # Normalize column names to lowercase for case-insensitivity
        peptide_df.columns = peptide_df.columns.str.lower()

        # Check if the required columns exist after normalization
        if 'start' not in peptide_df.columns or 'end' not in peptide_df.columns:
            return None

        output_peptide_list = []
        for _, row in peptide_df.iterrows():
            start = str(row['start']).strip()
            end = str(row['end']).strip()
            output_peptide_list.append(f"{start}-{end}")

        return output_peptide_list
    except Exception as e:
        return None



class FlagsParser:
    """
    A class to parse a flags file and validate its contents.
    """

    def __init__(self, file_path: str):
        self.file_path = file_path
        self.config = configparser.ConfigParser(allow_no_value=True)
        self.flags = {}

    def _read_flags(self) -> bool:
        """
        Reads the flags file, handling missing files.
        """
        if not os.path.exists(self.file_path) or not os.path.isfile(self.file_path):
            return False

        try:
            with open(self.file_path, 'r') as f:
                config_string = '[flags]\n' + f.read()
            self.config.read_string(config_string)
            self.flags = dict(self.config.items('flags'))
            return True
        except configparser.Error as e:
            return False

    def _validate_and_parse(self) -> Optional[Dict[str, Any]]:
        """
        Validates and parses the flags, exiting on error.
        """
        validation_rules = {
            'protein_name_states': {
                'type': list,
                'parser': self._parse_list_or_string,  # Must parse the string into a list of lists
                'validator': lambda x: (
                        isinstance(x, list) and
                        all(
                            isinstance(pair, list) and
                            len(pair) == 2 and
                            all(isinstance(s, str) for s in pair)
                            for pair in x
                        )
                )
            },
            'protein_sequence': {
                'type': (list, str),
                'parser': self._parse_list_or_string,
                'validator': self._validate_protein_sequence
            },
            'ph': {
                'type': float,
                'parser': float,
                'validator': lambda x: isinstance(x, float)
            },
            'd20_saturation': {
                'type': float,
                'parser': float,
                'validator': lambda x: 0.0 <= x <= 1.0
            },
            'temp': {
                'type': float,
                'parser': float,
                'validator': lambda x: isinstance(x, float)
            },
            'file_type': {
                'type': str,
                'parser': str,
                'validator': lambda x: x in ["DynamX", "HDXworkbench", "HDExaminer", "BioPharma", "Byos", "Custom"]
            }
        }

        parsed_flags = {}
        for key, rules in validation_rules.items():
            value = self.flags.get(key)
            if value is None or value.strip() == "":
                parsed_flags[key] = None
                continue

            try:
                # Strip comments and quotes before parsing
                value = value.split('#')[0].strip().strip('"')

                parsed_value = rules['parser'](value)
                if not rules['validator'](parsed_value):
                    raise ValueError(f"Validation failed for '{key}' with value '{value}'.")
                parsed_flags[key] = parsed_value
            except (ValueError, SyntaxError) as e:
                print(f"Error parsing flags file: Invalid value for '{key}'. Details: {e}")
                raise e

        return parsed_flags

    def _parse_list_or_string(self, value: str) -> Union[List[str], str]:
        if value.strip().startswith('[') and value.strip().endswith(']'):
            try:
                return ast.literal_eval(value)
            except (ValueError, SyntaxError):
                raise SyntaxError(f"Invalid list format: {value}")
        return value

    def _validate_protein_sequence(self, sequences: Union[List[str], str]) -> bool:
        allowed_chars = "ACDEFGHIKLMNPQRSTVWY"

        if isinstance(sequences, str):
            sequences = [sequences]

        for seq in sequences:
            if not all(c.upper() in allowed_chars for c in seq):
                return False
        return True

    def parse(self) -> Optional[Dict[str, Any]]:
        if not self._read_flags():
            return None

        return self._validate_and_parse()


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
            print(modification)
            if modification == "":
                modification = "N/A"

            replicate = int(row["replicate"]) - 1

            peptide_name = f"{row['start']}-{peptide}-{modification}-{row['end']}-{replicate}"

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
                hx_data.peptides[protein][state][peptide_name][exposure] = {
                    "uptake": float(row['centroid']) * float(row['charge']),
                    "start": int(row['start']),
                    "end": int(row['end']),
                    "sequence": str(row['peptide']),
                    "PTM": modification,
                    "time": float(exposure),
                    "mod": "A",
                    "envelope": "",
                    "replicate": int(row["replicate"]) - 1
                }
    return hx_data


def _parse_HDExaminer(file_path: str, flags_info: None) -> HxmsData:
    """
    Parses an HDExaminer CSV file into a HxmsData object.
    """
    hx_data = HxmsData()
    if not os.path.exists(file_path) or not os.path.isfile(file_path):
        return None

    # Step 1: Parse the header row to get time points and map them
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
            # Step 2: Read the rest of the file into a DataFrame
            hdx_df = pd.read_csv(file_path, skiprows=1)
        except Exception as e:
            print(f"Error reading file: {e}")
            return None

        # Step 3: Identify key columns and map time point data
        hdx_df.columns = hdx_df.columns.str.strip()

        # Check for required headers
        required_headers = ["State", "Protein", "Sequence", "Start", "End", "Charge"]
        for header in required_headers:
            if header not in hdx_df.columns:
                print(f"Missing required header: {header}")
                return None

        # Find all columns that contain '#D' for deuteration levels
        d_cols = [col for col in hdx_df.columns if (('#D' in col) and ('right' not in col))]

        # Step 4: Iterate through each row (peptide) and parse data for each time point
        for _, row in hdx_df.iterrows():
            used_rep = {}
            for item in time_points:
                used_rep[item]=0

            state = str(row['State']).strip()
            protein = str(row['Protein']).strip()
            # If protein is empty or 'nan', use a default name
            if not protein or protein == "nan":
                protein = "protein"

            peptide = str(row['Sequence']).strip()
            start = int(row['Start'])
            end = int(row['End'])
            modification = "N/A"

            # Add protein and state to the hx_data object if not present
            if protein not in hx_data.proteins:
                hx_data.proteins.append(protein)
                hx_data.state[protein] = []
                hx_data.peptides[protein] = {}
            if state not in hx_data.state[protein]:
                hx_data.state[protein].append(state)
                hx_data.peptides[protein][state] = {}

            # Create a unique key for the peptide

            # Step 5: Iterate through the time point columns

            for i, d_col in enumerate(d_cols):
                # Check if the time point data exists for this row
                if pd.notna(row[d_col]):
                    time_point_str = time_points[i]
                    time_point = float(time_point_str) if time_point_str != "inf" else float('inf')
                    uptake = float(row[d_col])
                    used_rep[time_point_str]+=1
                    peptide_name = f"{row['Start']}-{row['Sequence']}-{modification}-{row['End']}-0"
                    #peptide_name = f"{row['Start']}-{row['Sequence']}-{modification}-{row['End']}-"
                    #peptide_name = peptide_name+str(used_rep[time_point_str]-1)
                    if peptide_name not in hx_data.peptides[protein][state]:
                        hx_data.peptides[protein][state][peptide_name] = {}

                    if 0.0 not in hx_data.peptides[protein][state][peptide_name]:
                        hx_data.peptides[protein][state][peptide_name][0.0] = {
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
                        }

                    if time_point not in hx_data.peptides[protein][state][peptide_name]:
                        hx_data.peptides[protein][state][peptide_name][time_point] = {
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
                        }
                    else:

                        hx_data.peptides[protein][state][peptide_name][time_point] = {
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
                        }

    return hx_data


def _parse_custom(file_path: str, flags_info: None) -> HxmsData:
    """
    Parses a custom CSV file into a HxmsData object.
    """
    hx_data = HxmsData()
    if not os.path.exists(file_path) or not os.path.isfile(file_path):
        return None
    custom_headers = ["Protein","State","Start","Stop","Replicate","Modification","Bimodel_group","Exposure","Uptake","Envelope"]
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
        modification = str(row["Modification"]) if pd.notna(row["Modification"]) else "N/A"
        replicate = str(row["Replicate"]) if (pd.notna(row["Replicate"]) and row["Replicate"] != "") else "0"
        peptide_name = f"{row['Start']}-N/A-{modification}-{row['Stop']}-{replicate}"
        envelope = ""
        if row['Envelope'] != "":
            print(row['Envelope'])
            try:
                if ";" in row['Envelope']:
                    parts = row['Envelope'].replace("\n","").split(";")
                    if len(parts) > 1:
                        envelope = ",".join(parts)
                        print(envelope)
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

        hx_data.peptides[protein][state][peptide_name][str(row["Exposure"])] = {
            "uptake": float(row["Uptake"]),
            "start": int(row['Start']),
            "end": int(row['Stop']),
            "PTM": modification,
            "time": float(row["Exposure"]),
            "mod": str(row["Bimodel_group"]) if str(row["Bimodel_group"]) != "" else "A",
            "envelope": envelope,
            "replicate": replicate
        }
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
        peptide_name = f"{row['Start']}-{row['Sequence']}-{modification}-{row['End']}-0"
        if protein not in hx_data.proteins:
            hx_data.proteins.append(protein)
            hx_data.state[protein] = []
            hx_data.peptides[protein] = {}
        if state not in hx_data.state[protein]:
            hx_data.state[protein].append(state)
            hx_data.peptides[protein][state] = {}
        if peptide_name not in hx_data.peptides[protein][state]:
            hx_data.peptides[protein][state][peptide_name] = {}

        hx_data.peptides[protein][state][peptide_name][str(row["Exposure"])] = {
            "uptake": float(row["Uptake"]),
            "start": int(row['Start']),
            "end": int(row['End']),
            "sequence": str(row['Sequence']),
            "PTM": modification,
            "time": float(row["Exposure"]),
            "mod": "A",
            "envelope": "",
            "replicate": 0
        }
    return hx_data


def _parse_byos(file_path: str, flags_info: None) -> HxmsData:
    """
    Parses a Byos CSV file into a HxmsData object.
    """
    hx_data = HxmsData()
    if not os.path.exists(file_path) or not os.path.isfile(file_path):
        return None
    byos_headers = ['Condition', 'StartAA', 'EndAA', 'Sequence(unformatted)', 'Calc.M', 'Apex Time(Posit)',
                    'ExchangeTime', 'Replicate', '% deuteration']
    try:
        byos_df = pd.read_csv(file_path)
    except Exception as e:
        print(f"Error reading file: {e}")
        return None
    for header in byos_headers:
        if header not in byos_df.columns:
            return None
    for _, row in byos_df.iterrows():
        protein = "blank"
        state = str(row['Condition'])
        modification = "N/A"
        replicate = int(row["Replicate"]) - 1
        peptide_name = f"{row['StartAA']}-{row['Sequence(unformatted)']}-{modification}-{row['EndAA']}-0"
        if " - " in str(row["% deuteration"]) or pd.isna(row["% deuteration"]):
            continue
        if protein not in hx_data.proteins:
            hx_data.proteins.append(protein)
            hx_data.state[protein] = []
            hx_data.peptides[protein] = {}
        if state not in hx_data.state[protein]:
            hx_data.state[protein].append(state)
            hx_data.peptides[protein][state] = {}
        if peptide_name not in hx_data.peptides[protein][state]:
            hx_data.peptides[protein][state][peptide_name] = {}
        hx_data.peptides[protein][state][peptide_name][str(row["ExchangeTime"])] = {
            "uptake": float(row["% deuteration"]) * (int(row['EndAA']) - int(row['StartAA']) + 1),
            "start": int(row['StartAA']),
            "end": int(row['EndAA']),
            "sequence": str(row['Sequence(unformatted)']),
            "PTM": modification,
            "time": float(row["ExchangeTime"]),
            "mod": "A",
            "envelope": "",
            "replicate": replicate
        }
    return hx_data


def _parse_biopharma(file_path: str, flags_info: None) -> HxmsData:
    """
    Parses a BioPharma CSV file into a HxmsData object.
    """
    hx_data = HxmsData()
    if not os.path.exists(file_path) or not os.path.isfile(file_path):
        return None
    biopharam_headers = ["Protein", "Start", "End", "Sequence", "Modification", "Fragment", "MaxUptake", "MHP", "State",
                         "Exposure", "Center", "Center SD", "Uptake", "Uptake SD", "RT", "RT SD"]
    try:
        biopharma_df = pd.read_csv(file_path)
    except Exception as e:
        print(f"Error reading file: {e}")
        return None
    for header in biopharam_headers:
        if header not in biopharma_df.columns:
            return None
    for _, row in biopharma_df.iterrows():
        protein = str(row['Protein'])
        state = str(row['State'])
        modification = str(row["Modification"]) if pd.notna(row["Modification"]) else "N/A"
        peptide_name = f"{row['Start']}-{row['Sequence']}-{modification}-{row['End']}-0"

        if protein not in hx_data.proteins:
            hx_data.proteins.append(protein)
            hx_data.state[protein] = []
            hx_data.peptides[protein] = {}
        if state not in hx_data.state[protein]:
            hx_data.state[protein].append(state)
            hx_data.peptides[protein][state] = {}
        if peptide_name not in hx_data.peptides[protein][state]:
            hx_data.peptides[protein][state][peptide_name] = {}

        hx_data.peptides[protein][state][peptide_name][str(row["Exposure"])] = {
            "uptake": float(row["Uptake"]),
            "start": int(row['Start']),
            "end": int(row['End']),
            "sequence": str(row['Sequence']),
            "PTM": modification,
            "time": float(row["Exposure"]),
            "mod": "A",
            "envelope": "",
            "replicate": 0
        }
    return hx_data

def validate_and_parse_hxms_file(hxms_file_path:str)->[bool,str,HxmsData, dict]:
    if not os.path.exists(hxms_file_path):
        return [False,"File does not exist",None, None]
    if not os.path.isfile(hxms_file_path):
        return [False,"File is not a file",None, None]
    header_info = {}
    tp_info = {}
    ptm_info = {}
    skip_lines_TP = 0
    skip_lines_TP_done = False
    skip_lines_PTM = 0
    skip_lines_PTM_done = False
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

                        tp_info[data[2]+"-"+data[3]+"-"+data[4]+"-"+data[5]+"-"+data[7]] ={

                            "INDEX":int(data[1]),
                            "MOD": data[2],
                            "START": int(data[3]),
                            "END": int(data[4]),
                            "REP": int(data[5]),
                            "PTM_ID": int(data[6]),
                            "TIME(Sec)": float(data[7]),
                            "UPTAKE": float(data[8]),
                            "ENVELOPE": envelope,

                        }
                    except:
                        return [False, f"Incorrect TP info for {line}", None]

            if "TITLE_PTM" in line:
                if "TITLE_PTM  PTM_ID  CONTENT" not in line:
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
                        return [False, f"Incorrect TP info for {line}", None, None]

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

    for tp in tp_info:
        try:
            tp = tp_info[tp]
            peptide_name = f"{tp['START']}-N/A-{tp['PTM_ID']}-{tp['END']}"
            if peptide_name not in hx_data.peptides[protein][state]:
                hx_data.peptides[protein][state][peptide_name] = {}
            if str(tp['TIME(Sec)']) not in hx_data.peptides[protein][state][peptide_name]:
                hx_data.peptides[protein][state][peptide_name][str(tp['TIME(Sec)'])] = []
            hx_data.peptides[protein][state][peptide_name][str(tp['TIME(Sec)'])].append({
                "uptake": tp['UPTAKE'],
                "start": tp['START'],
                "end": tp['END'],
                "PTM": ptm_info[tp['PTM_ID']],
                "time": tp["TIME(Sec)"],
                "mod": tp["MOD"],
                "envelope": tp["ENVELOPE"],
                "replicate": tp["REP"]
            })
        except:
            return [False, "TP could not be read or converted properly", None]
    return [True,"",hx_data,ptm_info]
