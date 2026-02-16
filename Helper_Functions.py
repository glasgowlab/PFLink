from HXMS_IO import HxmsData
import csv
from typing import Optional
import zipfile
from pigeon_feather.spectra import *
from pigeon_feather.tools import custom_pad
from pigeon_feather.hxio import *
import re

def _create_peptide_chemical_formula(peptide:str)->dict:
    total_formula = {}
    residue_formulas = {
        'A': {'C': 3, 'H': 5, 'N': 1, 'O': 1},  # Alanine
        'R': {'C': 6, 'H': 12, 'N': 4, 'O': 1},  # Arginine
        'N': {'C': 4, 'H': 6, 'N': 2, 'O': 2},  # Asparagine
        'D': {'C': 4, 'H': 5, 'N': 1, 'O': 3},  # Aspartic Acid
        'C': {'C': 3, 'H': 5, 'N': 1, 'O': 1, 'S': 1},  # Cysteine
        'Q': {'C': 5, 'H': 8, 'N': 2, 'O': 2},  # Glutamine
        'E': {'C': 5, 'H': 7, 'N': 1, 'O': 3},  # Glutamic Acid
        'G': {'C': 2, 'H': 3, 'N': 1, 'O': 1},  # Glycine
        'H': {'C': 6, 'H': 7, 'N': 3, 'O': 1},  # Histidine
        'I': {'C': 6, 'H': 11, 'N': 1, 'O': 1},  # Isoleucine
        'L': {'C': 6, 'H': 11, 'N': 1, 'O': 1},  # Leucine
        'K': {'C': 6, 'H': 12, 'N': 2, 'O': 1},  # Lysine
        'M': {'C': 5, 'H': 9, 'N': 1, 'O': 1, 'S': 1},  # Methionine
        'F': {'C': 9, 'H': 9, 'N': 1, 'O': 1},  # Phenylalanine
        'P': {'C': 5, 'H': 7, 'N': 1, 'O': 1},  # Proline
        'S': {'C': 3, 'H': 5, 'N': 1, 'O': 2},  # Serine
        'T': {'C': 4, 'H': 7, 'N': 1, 'O': 2},  # Threonine
        'W': {'C': 11, 'H': 10, 'N': 2, 'O': 1},  # Tryptophan
        'Y': {'C': 9, 'H': 9, 'N': 1, 'O': 2},  # Tyrosine
        'V': {'C': 5, 'H': 9, 'N': 1, 'O': 1}  # Valine
    }
    water_formula = {'H': 2, 'O': 1}
    for AA in peptide.upper():
        if AA not in residue_formulas:
            continue
        residue = residue_formulas[AA]
        for element, count in residue.items():
            if element not in total_formula:
                total_formula[element] = 0
            total_formula[element] = total_formula[element] + count
    for element, count in water_formula.items():
        total_formula[element] = total_formula[element] + count
    return total_formula

def _convert_dformula_to_dict(dformula:str)->dict:
    pattern = re.compile(r'([A-Z][a-z]*)(\d*)')
    initial_counts = {}
    matches = pattern.findall(dformula)
    for element, count in matches:
        count = int(count) if count else 1
        initial_counts[element] = initial_counts.get(element, 0) + count
    final_formula = initial_counts.copy()
    if 'D' in final_formula:
        deuterium_count = final_formula.pop('D')
        final_formula['H'] = final_formula.get('H', 0) + deuterium_count

    return final_formula

def _subtract_chemical_formula_dicts(form1:dict,form2:dict)->str:
    diff = ""
    used = []
    for key in form1:
        if key in form2:
            used.append(key)
            difference = math.fabs(form1[key]-form2[key])
            if difference != 0:
                diff = diff + key + str(int(math.fabs(form1[key]-form2[key])))
        else:
            diff = diff + key + str(form1[key])
    for key in form2:
        if key not in used:
            diff = diff + key + str(form2[key])
    return diff


def _load_raw_ms_to_hdxms_data(hdxms_data, raw_spectra_path,protein_state:str = None):

    for state in hdxms_data.states:
        if protein_state != None:

            if state.state_name != protein_state:
                continue
        state_raw_spectra_path = os.path.join(raw_spectra_path, state.state_name)
        path_dict = {}
        folders = sorted(glob(state_raw_spectra_path + "/*"))

        for folder in folders:
            start, end, seq = os.path.basename(folder).split("-")
            start, end, seq = int(start), int(end), seq
            pep_idf = f"{start}-{end} {seq}"
            # pep_idf = f'{start+hdxms_data.n_fastamides}-{end}'
            path_dict[pep_idf] = folder

        for peptide in state.peptides:
            try:
                pep_sub_folder = path_dict[peptide.identifier]
            except:
                continue
            for tp in peptide.timepoints:
                if tp.deut_time == np.inf:
                    tp.isotope_envelope = None
                    continue
                elif tp.deut_time == 0:
                    try:
                        csv_file_path = glob(
                            f"{pep_sub_folder}/Non-D-*-z{tp.charge_state}.csv"
                        )[0]
                    except:
                        continue
                else:
                    csv_name = f"{int(tp.deut_time)}s-1-z{tp.charge_state}.csv"
                    csv_file_path = os.path.join(pep_sub_folder, csv_name)
                try:
                    df = tp.load_raw_ms_csv(csv_file_path, save_match=True)
                except:
                    continue

def _conver_PFhxms_to_hxms(dataset_file_path,protein_sequence,saturation,ph,temperature,envelope_file_path,protein_name,protein_state:str=None):
    hxms = HxmsData()
    hdxms_data_list = []
    cleaned = read_hdx_tables([dataset_file_path], [dataset_file_path], exclude=False, states_subset=[protein_state])
    hdxms_data = load_dataframe_to_hdxmsdata(
        cleaned,
        n_fastamides=0,
        protein_sequence=protein_sequence,
        fulld_approx=False,
        saturation=saturation,
        pH=ph,
        temperature=temperature
    )
    import zipfile
    with zipfile.ZipFile(envelope_file_path, 'r') as zip_ref:
        zip_ref.extractall(os.path.dirname(os.path.normpath(envelope_file_path)))
    raw_spectra_path = os.path.join(os.path.dirname(os.path.normpath(envelope_file_path)),
                                    os.path.splitext(os.path.basename(envelope_file_path))[0])

    _load_raw_ms_to_hdxms_data(
        hdxms_data,
        raw_spectra_path,
        protein_state
    )

    hdxms_data_list.append(hdxms_data)
    state_names = set([state.state_name for data in hdxms_data_list for state in data.states])
    if protein_name not in hxms.proteins:
        hxms.proteins.append(protein_name)
        hxms.state[protein_name] = []
        hxms.peptides[protein_name] = {}

    for state_name in state_names:
        protein_states = [state for data in hdxms_data_list for state in data.states if
                          state.state_name == state_name]
        all_peptides = [pep for state in protein_states for pep in state.peptides]
        all_timepoints = [tp for pep in all_peptides for tp in pep.timepoints]
        tp_dict = {}
        for tp_idx, tp in enumerate(all_timepoints):
            tp_dict[
                protein_name + "_" + state_name + "_" + str(tp.peptide.start) + "_" + str(tp.peptide.end) + "_" + str(
                    tp.deut_time)] = \
                tp.isotope_envelope is not None
        if state_name not in hxms.state[protein_name]:
            hxms.state[protein_name].append(state_name)
            hxms.peptides[protein_name][state_name] = {}

        for tp_idx, tp in enumerate(all_timepoints):
            peptide_name = f"{tp.peptide.start}~{tp.peptide.sequence}~N/A~{tp.peptide.end}~0"
            if peptide_name not in hxms.peptides[protein_name][state_name]:
                hxms.peptides[protein_name][state_name][peptide_name] = {}
            if tp.isotope_envelope is None and tp.deut_time != np.inf:
                continue
            key = protein_name + "_" + state_name + "_" + str(tp.peptide.start) + "_" + str(tp.peptide.end) + "_" + str(0.0)
            if key in tp_dict:
                if not tp_dict[key]:
                    continue
            else:
                continue
            if str(tp.deut_time) not in hxms.peptides[protein_name][state_name][peptide_name]:
                hxms.peptides[protein_name][state_name][peptide_name][str(tp.deut_time)] = [
                 {
                "uptake": float(tp.num_d),
                "start": int(tp.peptide.start),
                "end": int(tp.peptide.end),
                "sequence": str(tp.peptide.sequence),
                "PTM": "0000",
                "time": float(tp.deut_time),
                "mod": "A",
                "envelope": f"{','.join(f'{val:.3f}' for val in custom_pad(tp.isotope_envelope[:50].copy(), 50)) if tp.deut_time != np.inf else ''}",
                "replicate": 0,
                "match": f"{','.join(f'{val:.4f}' for val in tp.match) if tp.deut_time != np.inf else ''}",
                "mono_mz": f"{tp.mono_mz:.6f}" if tp.deut_time != np.inf else '',
                "charge_state": tp.charge_state,
                "RT": f"{tp.peptide.RT:.3f}",
                "match": f"{','.join(f'{val:.6f}' for val in tp.match) if tp.deut_time != np.inf else ''}",
                "raw_ms_fine_structure": ",".join(";".join(f"{m}:{i}" for m, i in zip(*val)) for val in tp.raw_ms_fine_structure) if tp.deut_time != np.inf else '',
                "score": f"{tp.score:.2f}"
            }]  
                # if tp.deut_time != np.inf:
                # print(tp.peptide.identifier, tp.deut_time)
                # print(tp.raw_ms_fine_structure)

            else:
                hxms.peptides[protein_name][state_name][peptide_name][str(tp.deut_time)].append(
                    {
                        "uptake": float(tp.num_d),
                        "start": int(tp.peptide.start),
                        "end": int(tp.peptide.end),
                        "sequence": str(tp.peptide.sequence),
                        "PTM": "0000",
                        "time": float(tp.deut_time),
                        "mod": "A",
                        "envelope": f"{','.join(f'{val:.3f}' for val in custom_pad(tp.isotope_envelope[:50].copy(), 50)) if tp.deut_time != np.inf else ''}",
                        "replicate": 0,
                        "match": f"{','.join(f'{val:.6f}' for val in tp.match) if tp.deut_time != np.inf else ''}",
                        "mono_mz": f"{tp.mono_mz:.6f}" if tp.deut_time != np.inf else '',
                        "charge_state": tp.charge_state,
                        "RT": f"{tp.peptide.RT:.3f}",
                        "match": f"{','.join(f'{val:.6f}' for val in tp.match) if tp.deut_time != np.inf else ''}",
                        "raw_ms_fine_structure": ",".join(";".join(f"{m}:{i}" for m, i in zip(*val)) for val in tp.raw_ms_fine_structure) if tp.deut_time != np.inf else '',
                        "score": f"{tp.score:.2f}"
                    })

    return hxms

def _check_peptide_list_format(file_path: str)->bool:
    if not os.path.exists(file_path) or not os.path.isfile(file_path):
        return False
    peptide_df = pd.read_csv(file_path)
    peptide_df.columns = peptide_df.columns.str.lower()
    if "start" not in peptide_df.columns:
        return False
    if "end" not in peptide_df.columns:
        return False
    return True



def _check_envelope_file(zip_file_path: str)->bool:
    if not os.path.exists(zip_file_path) or not os.path.isfile(zip_file_path):
        return False
    protein_states = []

    try:
        with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
            all_zip_entries = zip_ref.namelist()
            for entry in all_zip_entries:
                if entry.endswith('.csv'):
                    if ".csv" not in entry:
                        continue
                    folder_path = os.path.dirname(os.path.normpath(entry))
                    folder_name = os.path.basename(folder_path)
                    if len(folder_name.split("-")) != 3:
                        continue
                    state_path = os.path.dirname(os.path.normpath(folder_path))
                    state_name = os.path.basename(os.path.normpath(state_path))
                    if state_name not in protein_states:
                        protein_states.append(state_name)

    except zipfile.BadZipFile:
        print(f"Error: The file '{zip_file_path}' is not a valid zip file.")
        return False
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return False
    if len(protein_states) != 0:
        return True
    return False
def _sort_hxms_data_cc(hx_data: HxmsData) -> HxmsData:

    hx_data.proteins.sort()
    for protein in hx_data.state:
        hx_data.state[protein].sort()

    for protein in hx_data.peptides:
        for state in hx_data.peptides[protein]:
            peptide_data = hx_data.peptides[protein][state]
            for peptide_name, time_data in peptide_data.items():
                for time_key in time_data:
                    time_data[time_key].sort(
                        key=lambda measurement: (
                            measurement['replicate'],
                            measurement['time']
                        )
                    )
                peptide_data[peptide_name] = dict(sorted(
                    time_data.items(),
                    key=lambda item: float('inf') if item[0] == 'inf' else float(item[0])
                ))
            hx_data.peptides[protein][state] = dict(sorted(
                peptide_data.items(),
                key=lambda item: (
                    int(item[0].split('~')[0]),
                    int(item[0].split('~')[-1]),
                    item[0].split('~')[2],
                )
            ))

    return hx_data



def _get_avg_zero_uptake(hxms_data: HxmsData, protein: str, state: str, peptide_name: str) -> Optional[float]:

    zero_uptakes = []
    base_peptide_name = "~".join(peptide_name.split("~")[:-1])
    for p_name, p_data in hxms_data.peptides[protein][state].items():
        if "~".join(p_name.split("~")[:-1]) == base_peptide_name:
            if "0" in p_data:
                zero_uptakes.append(p_data["0"][0]["uptake"])
            if "0.0" in p_data:
                zero_uptakes.append(p_data["0.0"][0]["uptake"])
    if zero_uptakes:
        return sum(zero_uptakes) / len(zero_uptakes)
    return None


def check_file_format(file_path: str) -> str:

    known_formats = {
        "HDXworkbench": [
            "peptide", "charge", "start", "end", "experiment",
            "sample", "timepoint", "replicate", "discarded_replicate",
            "percentd_replicate", "centroid"
        ],
        "DynamX": [
            "Protein", "Start", "End", "Sequence", "Modification", "Fragment",
            "MaxUptake", "MHP", "State", "Exposure", "Center", "Center SD",
            "Uptake", "Uptake SD", "RT", "RT SD"
        ],
        "BioPharma": [
            "Retention Time", "Peptide", "Charge", "Residue From", "Residue To", "Ref", "0% Control",
            "100% Control", "stdev"
        ],
        "HDExaminer":[
            "State","Protein","Sequence","Search RT","Charge","Max D"
        ],
        "Custom": [
            "Protein","State","Start","Stop","Replicate","Modification","Bimodal_group","Exposure","Uptake","Envelope","RT","Conf","Z","Mono_M","M/Z"
        ]
    }

    if not os.path.exists(file_path):
        return "unknown"
    header_row_index = 0
    with open(file_path, 'r', newline='') as file:
        reader = csv.reader(file)
        for i, row in enumerate(reader):
            if any(h in row for h in ["peptide", "Protein", "Condition","State","Peptide"]):
                header_row_index = i
                break

    try:
        df = pd.read_csv(file_path, skiprows=header_row_index, nrows=0)
        file_headers = set(df.columns.str.strip())
    except Exception as e:
        return "unknown"
    for file_type, required_headers in known_formats.items():
        required_set = set(required_headers)
        if required_set.issubset(file_headers):
            return file_type
    return "unknown"


def validate_protein_sequence(sequence: str) -> bool:
    allowed_chars = "ACDEFGHIKLMNPQRSTVWY"

    for char in sequence:
        if not char.isalpha():
            return False
        if char.upper() not in allowed_chars:
            return False
    return True

def combine_hxms_data(hxms_data1:HxmsData,ptm_1_dic:dict,hxms_data2:HxmsData,ptm_2_dic:dict):
    largest_replicate = 0
    out_hxms = HxmsData()
    out_hxms.metadata = hxms_data1.metadata
    protein = hxms_data1.proteins[0]
    state = hxms_data1.state[protein][0]
    peptides = hxms_data1.peptides[protein][state]
    if protein not in out_hxms.proteins:
        out_hxms.proteins.append(protein)
        out_hxms.state[protein] = []
        out_hxms.peptides[protein] = {}
    if state not in out_hxms.state[protein]:
        out_hxms.state[protein].append(state)
        out_hxms.peptides[protein][state] = {}


    for peptide in peptides:
        if peptide not in out_hxms.peptides[protein][state]:
            out_hxms.peptides[protein][state][peptide] = {}
        for tp in peptides[peptide]:
            if tp not in out_hxms.peptides[protein][state][peptide]:
                out_hxms.peptides[protein][state][peptide][tp] = []
            for measurment in peptides[peptide][tp]:
                try:
                    rep_value = int(measurment["replicate"])
                    if rep_value > largest_replicate:
                        largest_replicate = rep_value
                    out_hxms.peptides[protein][state][peptide][tp].append(measurment)
                except:
                    return
    largest_replicate += 1
    peptides = hxms_data2.peptides[protein][state]
    for peptide in peptides:
        if peptide not in out_hxms.peptides[protein][state]:
            out_hxms.peptides[protein][state][peptide] = {}
        for tp in peptides[peptide]:
            if tp not in out_hxms.peptides[protein][state][peptide]:
                out_hxms.peptides[protein][state][peptide][tp] = []
            for measurment in peptides[peptide][tp]:
                try:
                    measurment["replicate"] = measurment["replicate"] + largest_replicate
                    out_hxms.peptides[protein][state][peptide][tp].append(measurment)
                except:
                    return

    out_hxms = _sort_hxms_data_cc(out_hxms)
    return out_hxms
