import gradio as gr
import pandas as pd
import tempfile
import os


# This function contains the core logic for processing the Centroid data.
# It is called by process_centroid_data after all inputs are validated.
def _run_data_processing(
        temperature: float,
        ph: float,
        saturation: float,
        file_format: str,
        dataset_file: gr.File,
        protein_name: str,
        protein_state_name: str,
        protein_sequence: str,
        protein_state_rename:str,
        protein_name_rename:str,
        peptide_list_parsed: list,
        include_exclude_bool: bool,
        envelope_file:  gr.File,
        replicate_ID: int,
        bimodel_group: str
) -> (str, str):
    """
    This function processes the user's inputs and files for the Centroid tab.
    It assumes all inputs have been validated and now creates a temporary
    output file for the user to download.
    """
    parser_map = {
        "DynamX": _parse_dynamX,
        "HDXworkbench": _parse_HDXWorkbench,
        "BioPharma": _parse_biopharma,
        "Byos": _parse_byos,
        "HDExaminer": _parse_HDExaminer,
        "Custom": _parse_custom,
    }

    # Create a fake flags dictionary from the user inputs
    flags = {
        'protein_name_states': [[protein_name,protein_state_name]],
        'protein_sequence': protein_sequence,
        'ph': ph,
        'd20_saturation': saturation,
        'temp': temperature,
        'file_type': file_format
    }

    dataset_file_path = dataset_file.name

    if envelope_file:
        hxms = _conver_PFhxms_to_hxms(dataset_file_path,protein_sequence,saturation,ph,temperature,envelope_file.name,protein_name)
    else:
        hxms = parser_map[file_format](dataset_file_path, flags)

    # Sanitize names for use in filename
    sanitized_protein_name = protein_name.replace(" ", "_")
    sanitized_state_name = protein_state_name.replace(" ", "_")

    # Create a temporary file with a descriptive name to store the HXMS output
    temp_output_file = tempfile.NamedTemporaryFile(
        prefix=f"{sanitized_protein_name}_{sanitized_state_name}_",
        suffix=".hxms",
        delete=False
    )
    temp_output_path = temp_output_file.name
    temp_output_file.close()

    # Write the data to the temporary file
    if hxms:
        hxms_good = write_hxms_file(hxms, temp_output_path, flags,protein_state_rename, protein_name_rename,peptide_list_parsed,include_exclude_bool,replicate_ID,bimodel_group)
        if not hxms_good:
            raise gr.Error(
                "Failed to parse your input data. "
                "Please make sure the protein name and state match your input data or ensure your input data is valid!"
                "Also, check your peptide inclusion/exclusion list.",
                print_exception=False)

        else:
            gr.Info("Your hxms file is ready to download at the bottom!")

    format_found = check_file_format(dataset_file.name)

    try:
        file_size_bytes = os.path.getsize(temp_output_path)
        if file_size_bytes/1028 < 1:
            gr.Warning("Your results file is rather small. Please check your inputs for errors.",duration=30)
    except:
        gr.Warning("Your results file is rather small. Please check your inputs for errors.",duration=30)

    confirmation_message = (
        f"Settings:\n"
        f"  - Temperature: {temperature}\n"
        f"  - pH: {ph}\n"
        f"  - Saturation: {saturation}\n"
        f"  - File Format: {file_format}\n\n"
        f"Protein and State Data:\n"
        f"  - Name: {protein_name}\n"
        f"  - State: {protein_state_name}\n"
        f"  - Sequence: {protein_sequence}\n\n"
        f"Files:\n"
        f"  - Dataset File: {'Uploaded' if dataset_file else 'Not Uploaded'}\n"
    )
    confirmation_message+=f"Centroid Data Received Successfully!\n\n"
    if (envelope_file is not None) and (envelope_file.name != ""):
        confirmation_message += f"Envelope Data Received Successfully!\n\n"
    if format_found != file_format:
        confirmation_message+=f"The file may have an error!\n"
        confirmation_message += f"The detected file format was {format_found} while you selected {file_format}\n"
    # Return both the confirmation message and the path to the downloadable file
    return confirmation_message, temp_output_path


# Define the function for the Centroid tab processing with input validation.
def process_centroid_data(
        temperature: str,
        ph: str,
        saturation: str,
        file_format: str,
        dataset_file: gr.File,
        protein_name: str,
        protein_state_name: str,
        protein_sequence: str,
        protein_names: str,
        protein_state_names: str,
        use_text_input:bool,
        protein_state_rename:str,
        protein_name_rename:str,
        include_exclude: str,
        peptide_list: gr.File,
        envelope_file: gr.File,
        replicate_ID: str,
        bimodel_group: str

) -> (str, str):
    """
    This function first validates all required inputs for the Centroid tab.
    If all inputs are valid, it calls the core processing function.
    """
    if not use_text_input:
        protein_state_name = protein_state_names
        protein_name = protein_names

    # 1. Check for required file uploads
    if not dataset_file:
        gr.Warning("Please upload a Dataset File.",duration=30)
        return "Process failed: Dataset File is missing.", None

    if envelope_file:
        if not _check_envelope_file(envelope_file.name):
            gr.Warning("Bad envelope file. Please check the example file format!",duration=30)
            return "Process failed: Bad envelope file.", None
        if file_format != "HDExaminer":
            gr.Warning("Envelope data is only allowed for HDExaminer at this time.",duration=30)
            return "Process failed: Non HDExaminer envelope data.", None

    # 2. Check for required text inputs
    if not temperature:
        gr.Warning("Please provide a value for Temperature.",duration=30)
        return "Process failed: Temperature is missing.", None
    if not ph:
        gr.Warning("Please provide a value for pH.",duration=30)
        return "Process failed: pH is missing.", None
    if not saturation:
        gr.Warning("Please provide a value for Saturation.",duration=30)
        return "Process failed: Saturation is missing.", None
    if not protein_name:
        gr.Warning("Please provide a value for Protein Name.",duration=30)
        return "Process failed: Protein Name is missing.", None
    if not protein_state_name:
        gr.Warning("Please provide a value for Protein State Name.",duration=30)
        return "Process failed: Protein State Name is missing.", None
    if not protein_sequence:
        gr.Warning("Please provide a value for Protein Sequence.",duration=30)
        return "Process failed: Protein Sequence is missing.", None

    # 3. Validate numeric and sequence inputs
    try:
        temp_float = float(temperature)
        ph_float = float(ph)
        saturation_float = float(saturation)
    except ValueError:
        gr.Error("Temperature, pH, and Saturation must be numeric.")
        return "Process failed: Non-numeric input for experimental settings.", None

    if not (200 <= temp_float <= 370):
        gr.Warning("Temperature must be between 200 and 370 Kelvin.",duration=30)
        return "Process failed: Invalid Temperature value.", None
    if not (0 <= ph_float <= 14):
        gr.Warning("pH must be between 0 and 14.",duration=30)
        return "Process failed: Invalid pH value.", None
    if not (0 <= saturation_float <= 1):
        gr.Warning("Saturation must be between 0 and 1.",duration=30)
        return "Process failed: Invalid Saturation value.", None

    if not validate_protein_sequence(protein_sequence):
        gr.Warning("Invalid protein sequence. Please use the valid 20 AA chars.",duration=30)
        return "Process failed: Invalid protein sequence.", None

    peptide_list_parsed = None
    if peptide_list is not None:
        if not _check_peptide_list_format(peptide_list.name):
            gr.Warning("Invalid peptide list! Please download the example csv for guidance on the format!",duration=30)
            return "Process failed: Invalid peptide list.", None
        else:
            peptide_list_parsed = _parse_peptide_list(peptide_list.name)

    if (replicate_ID is not None) and (replicate_ID != ""):
        try:
            replicate_ID = int(replicate_ID)
            if replicate_ID < 0:
                gr.Warning("Replicate ID must be an integer bigger or equal to 0.",duration=30)
                return "Process failed: Replicate ID is not a int bigger or equal to 0.", None
        except:
            gr.Warning("Replicate ID must be an integer.",duration=30)
            return "Process failed: Replicate ID is not integer.",None
    else:
        replicate_ID = None

    if (bimodel_group is not None) and (bimodel_group != ""):
        if len(bimodel_group) != 1:
            gr.Warning("Bimodel Group must be a single capital letters only.",duration=30)
            return "Process failed: Bimodel Group is not a single letter only.",None
        if not bimodel_group.isalpha():
            gr.Warning("Bimodel Group must be capital letters only.",duration=30)
            return "Process failed: Bimodel Group is not a single letter only.",None
        if bimodel_group.upper() != bimodel_group:
            gr.Warning("Bimodel Group must be capital letters only.",duration=30)
            return "Process failed: Bimodel Group is not a single capital letter only.",None
    else:
        bimodel_group = None

    if include_exclude == "Include":
        include_exclude_bool = True
    else:
        include_exclude_bool = False

    format_found = check_file_format(dataset_file.name)

    if format_found != file_format:
        gr.Warning("The file format you selected and the one you upload may be different formats! Please check.",duration=30)

    # 4. If all validations pass, call the core processing function
    return _run_data_processing(
        temp_float,
        ph_float,
        saturation_float,
        file_format,
        dataset_file,
        protein_name,
        protein_state_name,
        protein_sequence,
        protein_state_rename,
        protein_name_rename,
        peptide_list_parsed,
        include_exclude_bool,
        envelope_file,
        replicate_ID,
        bimodel_group
    )


#


# This function provides real-time validation feedback using gr.Warning and gr.Error.
def validate_inputs_experimental(temp: str, ph: str, saturation: str):
    try:
        if (temp is not None) and (temp != ""):
            float(temp)
        if (ph is not None) and (ph != ""):
            float(ph)
        if (saturation is not None) and (saturation != ""):
            float(saturation)
    except:
        gr.Warning("Temperature, ph, and saturation must be numeric",duration=30)
        return

    if (temp is not None) and (temp != ""):
        temp = float(temp)
        if not (200 <= temp <= 370):
            gr.Warning("Temperature must be between 200 and 370 Kelvin.",duration=30)
            return

    if (ph is not None) and (ph != ""):
        ph = float(ph)
        if not (0 <= ph <= 14):
            gr.Warning("pH must be between 0 and 14.",duration=30)
            return

    if (saturation is not None) and (saturation != ""):
        saturation = float(saturation)
        if not (0 <= saturation <= 1):
            gr.Warning("Saturation must be between 0 and 1.",duration=30)
            return

    if temp is None or ph is None or saturation is None:
        gr.Warning("Please give a numeric value to temp, ph, and saturation.",duration=30)
        return


def validate_inputs_protein(protein_sequence: str):
    if (protein_sequence is not None) and protein_sequence != "":
        if not validate_protein_sequence(protein_sequence):
            gr.Warning("Invalid protein sequence. Please use the valid 20 AA chars.",duration=30)
            return


def process_peptide_list(peptide_file: gr.File):
    if peptide_file is None:
        return None
    if not _check_peptide_list_format(peptide_file.name):
        gr.Warning("Invalid peptide list! Please download the example csv for guidance on the format!",duration=30)
def envelope_file_file_check(envelope_file: gr.File) -> None:
    if envelope_file is None:
        return None
    if not _check_envelope_file(envelope_file.name):
        gr.Warning("Invalid Envelope file! Please see the example file format.",duration=30)
# Added function to process flags file
def process_flags(flags_file: gr.File):
    """
    This function processes the flags file immediately upon upload.
    It returns a message and updates the text boxes with the parsed values.
    """
    if flags_file is None:
        return (
            "No flags file uploaded.",
            gr.update(value=None),
            gr.update(value=None),
            gr.update(value=None),
            gr.update(value=None),
            gr.update(value=None),
            gr.update(value=None),
            gr.update(value=None),
            gr.update(visible=False),
            gr.update(visible=False),
            gr.update(value=True)
        )
    else:
        try:
            parser = FlagsParser(flags_file.name)
            flags = parser.parse()

            # The flags file may contain a list of proteins and states, so we'll just use the first one for the text boxes.
            protein_state_info = flags.get('protein_name_states', [None])[0]
            if protein_state_info:
                # Split the string to get name and state, handling cases where it might be a single string
                parts = protein_state_info[0]
                protein_name_val = parts[0]
                protein_state_name_val = parts[1]
            else:
                protein_name_val = None
                protein_state_name_val = None

            protein_sequence_val = flags.get('protein_sequence', [None])[0]
            temp_val = flags.get('temp')
            ph_val = flags.get('ph')
            sat_val = flags.get('d20_saturation')
            file_format_val = flags.get('file_type')

            return (
                "Flags file processed successfully!",
                gr.update(value=protein_name_val, visible=True),
                gr.update(value=protein_state_name_val, visible=True),
                gr.update(value=protein_sequence_val),
                gr.update(value=temp_val),
                gr.update(value=ph_val),
                gr.update(value=sat_val),
                gr.update(value=file_format_val),
                gr.update(visible=False),  # Hide protein names dropdown
                gr.update(visible=False),  # Hide protein state names dropdown
                gr.update(value=True)
            )
        except Exception as e:
            gr.Warning(f"Error processing flags file: {e}")
            return (
                "Error processing flags file. Please see the popup for details.",
                gr.update(value=None),
                gr.update(value=None),
                gr.update(value=None),
                gr.update(value=None),
                gr.update(value=None),
                gr.update(value=None),
                gr.update(value=None),
                gr.update(visible=False),  # Hide protein names dropdown
                gr.update(visible=False),  # Hide protein state names dropdown
                gr.update(value=False)
            )


# This single function handles all logic for showing/hiding protein fields
def update_protein_fields(
        dataset_file: gr.File,
        use_text_input: bool,
        protein_name_text: str,
        protein_state_text: str
):
    """
    Manages the visibility and content of the protein textboxes and dropdowns.
    It prioritizes dropdowns if a file is present and the textboxes are empty.
    """
    file_uploaded = dataset_file is not None
    show_dropdowns = file_uploaded and not use_text_input
    show_textboxes = not show_dropdowns

    dropdown_names = []
    dropdown_states = []
    hxms = None
    format_found = "Dynamx"
    if file_uploaded:
        format_found = check_file_format(dataset_file.name)
        if format_found != "unknown":
            parser_map = {
                "DynamX": _parse_dynamX,
                "HDXworkbench": _parse_HDXWorkbench,
                "BioPharma": _parse_biopharma,
                "Byos": _parse_byos,
                "HDExaminer": _parse_HDExaminer,
                "Custom": _parse_custom,
            }
            hxms = parser_map[format_found](dataset_file.name, None)
            dropdown_names = hxms.proteins

            if dropdown_names:
                dropdown_states = hxms.state.get(dropdown_names[0], [])
        else:
            show_dropdowns = False
            show_textboxes = True

    return (
        # protein_name textbox
        gr.update(visible=show_textboxes),
        # protein_state_name textbox
        gr.update(visible=show_textboxes),
        # protein_names dropdown
        gr.update(visible=show_dropdowns, choices=dropdown_names, value=dropdown_names[0] if dropdown_names else None),
        # protein_state_names dropdown
        gr.update(visible=show_dropdowns, choices=dropdown_states,
                  value=dropdown_states[0] if dropdown_states else None),
        # Checkbox visibility and value
        gr.update(visible=file_uploaded, value=use_text_input),
        # Return the hxms object to be stored in the state
        hxms,
        gr.update(value=format_found)
    )


def update_states_dropdown(selected_protein, hxms_data):
    """
    Updates the protein state dropdown based on the selected protein name.
    """
    if hxms_data and selected_protein:
        states = hxms_data.state.get(selected_protein, [])
        return (
            gr.update(choices=states, value=states[0] if states else None)
        )
    return gr.update(choices=[], value=None)


def load_centroid_data_example():
    #dataset_file,peptide_list,flags_file,include_exclude
    return  (gr.update(value="./test_data/ecDHFR_tutorial.csv"),
             gr.update(value="./test_data/rangeslist.csv"),
             gr.update(value="Include"))

def load_enevlope_data_example():
    #dataset_file,envlope_data,peptide_list,flags_file,include_exclude
    return (gr.update(value="./test_data/ecDHFR_tutorial.csv"),
            gr.update(value="./test_data/SpecExport.zip"),
             gr.update(value="./test_data/rangeslist.csv"),
             gr.update(value="Include"))


def check_replicate_id(replicate_id: str):
    if (replicate_id is not None) and replicate_id != "":
        try:
            replicate_id = int(replicate_id)
            if replicate_id < 0:
                gr.Warning("Replicate ID must be an integer bigger or equal to 0.",duration=30)
                return
        except:
            gr.Warning("Replicate ID must be an integer.",duration=30)
            return
    return

def check_bimodel_group(bimodel_group: str):
    if (bimodel_group is not None) and bimodel_group != "":
        try:
            if len(bimodel_group) != 1:
                gr.Warning("Bimodel Group must a single capital letters only.",duration=30)
                return
            if not bimodel_group.isalpha():
                gr.Warning("Bimodel Group must be capital letters only.",duration=30)
                return
            if bimodel_group.upper() != bimodel_group:
                gr.Warning("Bimodel Group must be capital letters only.",duration=30)
                return
        except:
            gr.Warning("Bimodel Group must be capital letters only.",duration=30)
            return
    return

def delayed_flags_processing(flags_file):
    return gr.update(value="./test_data/flag-DHFR.txt")


def validate_HXMS(hxms_file):
    if (hxms_file is not None) and (hxms_file.name != ""):
        valid = validate_and_parse_hxms_file(hxms_file.name)
        if not valid[0]:
            gr.Warning("Invalid HXMS file. Please check and upload again", duration=30)
            gr.Warning(valid[1], duration=30)



def run_hxms_combine(hxms_1:gr.File,hxms_2:gr.File):
    if (hxms_1 is None) or (hxms_2 is None):
        gr.Warning("Please ensure you have uploaded both hxms files.", duration=30)
        return
    hxms_data_1 = validate_and_parse_hxms_file(hxms_1.name)
    hxms_data_2 = validate_and_parse_hxms_file(hxms_2.name)
    if (not hxms_data_1[0]) or (not hxms_data_2[0]):
        gr.Warning("One of your hxms files are wrong. Please fix them before trying to merge.", duration=30)
        return

    hxms_1_hxms_object: HxmsData = hxms_data_1[2]
    hxms_2_hxms_object: HxmsData = hxms_data_2[2]

    if hxms_1_hxms_object.proteins != hxms_2_hxms_object.proteins:
        gr.Warning("Your protein names do no match. Please sure you are using the same protein and try again.", duration=30)
        return
    if hxms_1_hxms_object.state != hxms_2_hxms_object.state:
        gr.Warning("Your protein states do no match. Please sure you are using the same state and try again.", duration=30)
        return
    if hxms_1_hxms_object.metadata['protein_sequence'] != hxms_2_hxms_object.metadata['protein_sequence'] :
        gr.Warning("Your protein sequences do no match. Please sure you are using the same sequence and try again.", duration=30)
        return
    if hxms_1_hxms_object.metadata['temperature'] != hxms_2_hxms_object.metadata['temperature'] :
        gr.Warning("Your temperatures do no match. Please sure you are using the same temperature and try again.", duration=30)
        return
    if hxms_1_hxms_object.metadata['ph'] != hxms_2_hxms_object.metadata['ph'] :
        gr.Warning("Your phs do no match. Please sure you are using the same phs and try again.", duration=30)
        return
    if hxms_1_hxms_object.metadata['saturation'] != hxms_2_hxms_object.metadata['saturation'] :
        gr.Warning("Your saturations do no match. Please sure you are using the same saturations and try again.", duration=30)
        return

    out_hxms = combine_hxms_data(hxms_1_hxms_object, hxms_data_1[-1], hxms_2_hxms_object, hxms_data_2[-1])

    from HXMS_IO import write_hxms_file_combined_test
    temp_output_file = tempfile.NamedTemporaryFile(
        prefix=f"{hxms_1_hxms_object.metadata['protein_name']}_{hxms_1_hxms_object.metadata['protein_state']}_combined_",
        suffix=".hxms",
        delete=False
    )
    temp_output_path = temp_output_file.name
    temp_output_file.close()

    print(temp_output_path)
    write_hxms_file_combined_test(out_hxms,temp_output_path)


    if temp_output_file is not None:
        return temp_output_path
    else:
        gr.Warning("Error, could not make HXMS file")
        return None



# Create the Gradio Interface with Tabs
with gr.Blocks(title="PFLink") as demo:
    gr.Markdown("# PFLink")
    gr.Markdown("PFLink takes a variety of HDX-MS file formats and converts it to the HXMS file format.")
    gr.Markdown(
        "HXMS files can be processed in [PFNet](https://pfnet-python.readthedocs.io/en/latest/) to predicts ΔGop for protein residues.")
    with gr.Tabs():
        with gr.TabItem("Create HXMS files"):
            gr.Markdown("---")
            gr.Markdown("**Basic instructions**")
            gr.Markdown(
                "Use the *Centroid* tab to upload your centroid data.\n"
                "Use the *Envelope* tab in addition to centroid data if you have envelope data (only with HDExaminer at this time).\n\n"
                "If you are using the 'Custom' data format, please upload it to the centroid data section, regardless of whether you have envelopes in the file.\n\n"
                "Each file upload only accepts one file. The online tool only processes one state/protein at a time within that file.\n\n"
                "If you have multiple state/protein(s) in that file, you must select the one you want to process.\n\n"
                "You can download example files using the 'Download example data' button below.\n\n"
                "The flags-file and custom-format CSV can be found inside the example data download zip.\n\n"
                "If the option does not say optional, then it is **required**.\n\n"
                "---")


            gr.Markdown("Example datasets")
            with gr.Row():
                submit_centroid_test_btn = gr.Button("Centroid data load example")
                submit_envelope_test_btn = gr.Button("Envelope data load example")
                download_data = gr.DownloadButton("Download example data", visible=True,value="./test_data/test_data.zip")
            gr.Markdown("---")
            data_state = gr.State()
            with gr.Tabs():
                with gr.TabItem("Centroid"):
                    gr.Markdown("### Centroid File Uploads (.csv)")
                    with gr.Column(elem_id="centroid-files-column"):
                        dataset_file = gr.File(label="Upload Dataset File", file_types=[".csv"])
                        # New output box to show the status of the flags file


                with gr.TabItem("Envelope"):
                    gr.Markdown("### Envelope File Upload (.zip)")
                    with gr.Column(elem_id="envelope-files-column"):
                        envelope_file = gr.File(label="Upload Envelope File", file_types=[".zip"])


            gr.Markdown("**Global settings**")
            gr.Markdown(
                "Flags file (.txt) - optional. A flags file can be used to autofill the global settings. It is not necessary, but it can speed up entering data.")
            flags_file = gr.File(label="Upload Flags File", file_types=[".txt"], height=120)
            flags_message = gr.Textbox(label="Flags File Status", lines=2, interactive=False)
            gr.Markdown("---")
            # Wrap these sections in columns to make the layout vertical
            with gr.Row():
                with gr.Column(scale=1):
                    gr.Markdown("Experimental settings")
                    temp = gr.Textbox(label="Temperature in Kelvin (float 200-350K)", placeholder="e.g., 278.15", scale=1,
                                      interactive=True)
                    ph = gr.Textbox(label="pH - (pH read) (0-14)", placeholder="e.g., 7.4", scale=1, interactive=True)
                    saturation = gr.Textbox(label="Saturation (0.0-1.0)", placeholder="e.g., 0.7", scale=1, interactive=True)
                with gr.Column(scale=1):
                    gr.Markdown("Protein information")
                    use_text_input = gr.Checkbox(label="Manually type protein name/state", value=False, visible=True)
                    # Dropdowns are hidden by default
                    protein_names = gr.Dropdown([], label="Protein names", scale=1, visible=False)
                    protein_state_names = gr.Dropdown([], label="Protein state names", scale=1, visible=False)
                    # Textboxes are visible by default
                    protein_name = gr.Textbox(label="Protein name (Case sensitive)", placeholder="e.g., ecDHFR", scale=1,
                                              interactive=True)
                    protein_state_name = gr.Textbox(label="Protein state name (Case sensitive)", placeholder="e.g., APO", scale=1,
                                                    interactive=True)
                    protein_sequence = gr.Textbox(label="Protein sequence (20 common AAs 1 letter codes only)", placeholder="e.g., MGGATS...", scale=1,
                                                  interactive=True)
            gr.Markdown("File settings - choose your HDX input file format.")
            with gr.Row():
                file_format = gr.Dropdown(["DynamX", "HDXworkbench", "HDExaminer", "BioPharma", "Byos", "Custom"],
                                          label="Formats", value="DynamX", scale=1)

            gr.Markdown("Protein rename and more - optional. Renames your state/protein in your output HXMS file. As well as forcing a replicate ID or bimodal group")
            with gr.Row():

                protein_state_rename = gr.Textbox(label="Protein state rename", placeholder="e.g., APO", scale=1,
                                                interactive=True)
                protein_name_rename = gr.Textbox(label="Protein name rename", placeholder="e.g., LacI", scale=1,
                                              interactive=True)
                replicate_ID = gr.Textbox(label="replicate ID (Must be positive integers only)", placeholder="e.g., 0", scale=1,
                                              interactive=True)
                bimodel_group = gr.Textbox(label="Bimodel group (Must be capital letters only)", placeholder="e.g., A", scale=1,
                                              interactive=True)

            gr.Markdown("Peptide inclusion/exclusion (.csv) - optional. Filters peptides that you want to include/exclude in your output HXMS file.")
            with gr.Row():
                include_exclude = gr.Dropdown(["Include", "Exclude"],label="Filter list type", value="Include", scale=1,interactive=True)
                peptide_list = gr.File(label="Upload Dataset File", file_types=[".csv"])



            gr.Markdown("---")
            gr.Markdown("Submit and output data")
            with gr.Column():
                submit_centroid_btn = gr.Button("Process Data")
                output_centroid_message = gr.Textbox(label="Result", lines=10)
                download_file = gr.File(label="Download HXMS File", interactive=False)

        with gr.TabItem("Combine HXMS files"):
            gr.Markdown("###HXMS file uploads###")
            with gr.Column(elem_id="hxms-files-column"):
                HXMS_file1 = gr.File(label="Upload the first HXMS file", file_types=[".hxms"])
                HXMS_file2 = gr.File(label="Upload the second HXMS file", file_types=[".hxms"])
                submit_HXMS_btn = gr.Button("Process Data")
                download_hxms_combine_file = gr.File(label="Download HXMS File", interactive=False)


    # Link the change events of the number inputs to the new validation function
    temp.change(fn=validate_inputs_experimental, inputs=[temp, ph, saturation], outputs=None)
    ph.change(fn=validate_inputs_experimental, inputs=[temp, ph, saturation], outputs=None)
    saturation.change(fn=validate_inputs_experimental, inputs=[temp, ph, saturation], outputs=None)

    protein_sequence.submit(fn=validate_inputs_protein, inputs=[protein_sequence], outputs=None)
    protein_sequence.change(fn=validate_inputs_protein, inputs=[protein_sequence], outputs=None)

    replicate_ID.submit(fn=check_replicate_id, inputs=[replicate_ID], outputs=None)
    replicate_ID.change(fn=check_replicate_id, inputs=[replicate_ID], outputs=None)

    bimodel_group.submit(fn=check_bimodel_group, inputs=[bimodel_group], outputs=None)
    bimodel_group.change(fn=check_bimodel_group, inputs=[bimodel_group], outputs=None)

    # New event handler to process the flags file and update the text boxes immediately
    flags_file.change(
        fn=process_flags,
        inputs=flags_file,
        outputs=[flags_message, protein_name, protein_state_name, protein_sequence, temp, ph, saturation, file_format,
                 protein_names, protein_state_names,use_text_input]
    )

    peptide_list.change(
        fn=process_peptide_list,
        inputs=peptide_list,
        outputs=None
    )

    dataset_file.change(
        fn=update_protein_fields,
        inputs=[dataset_file, use_text_input, protein_name, protein_state_name],
        outputs=[protein_name, protein_state_name, protein_names, protein_state_names, use_text_input, data_state,file_format]
    )

    # Event for the new checkbox to control the fields
    use_text_input.change(
        fn=update_protein_fields,
        inputs=[dataset_file, use_text_input, protein_name, protein_state_name],
        outputs=[protein_name, protein_state_name, protein_names, protein_state_names, use_text_input, data_state,file_format]
    )

    # New event listener to update the state dropdowns when a protein name is selected
    protein_names.change(
        fn=update_states_dropdown,
        inputs=[protein_names, data_state],
        outputs=[protein_state_names]
    )


    # Link the button click events to their respective functions
    submit_centroid_btn.click(
        fn=process_centroid_data,
        inputs=[
            temp, ph, saturation, file_format, dataset_file, protein_name, protein_state_name,
            protein_sequence, protein_names,protein_state_names,use_text_input,protein_state_rename,
            protein_name_rename,include_exclude,peptide_list,envelope_file,replicate_ID,bimodel_group],
        outputs=[output_centroid_message, download_file]
    )

    HXMS_file1.change(
        fn=validate_HXMS,
        inputs=[HXMS_file1],
        outputs=None
    )

    HXMS_file2.change(
        fn=validate_HXMS,
        inputs=[HXMS_file2],
        outputs=None
    )

    submit_HXMS_btn.click(
        fn=run_hxms_combine,
        inputs=[HXMS_file1,HXMS_file2],
        outputs=[download_hxms_combine_file]

    )

    envelope_file.change(
        fn=envelope_file_file_check,
        inputs=[envelope_file],
        outputs=None
    )


    submit_centroid_test_btn.click(
        fn=load_centroid_data_example,
        inputs=None,
        outputs=[dataset_file,peptide_list,include_exclude]
    ).then(
    fn=delayed_flags_processing,
    inputs=flags_file,
    outputs=[flags_file])




    submit_envelope_test_btn.click(
        fn=load_enevlope_data_example,
        inputs=None,
        outputs=[dataset_file,envelope_file,peptide_list,include_exclude]
    ).then(
    fn=delayed_flags_processing,
    inputs=flags_file,
    outputs=[flags_file])


# Launch the Gradio application
if __name__ == "__main__":
    try:
        from pigeon_feather.hxio import *
    except ImportError:
        os.system("pip install ./PIGEON-FEATHER/.")
        os.system("pip install --upgrade ./PIGEON-FEATHER/.")
    from HXMS_IO import write_hxms_file, HxmsData
    from Helper_Functions import validate_protein_sequence, check_file_format, _check_peptide_list_format, \
        _check_envelope_file, _conver_PFhxms_to_hxms,combine_hxms_data
    from Parsers import FlagsParser, _parse_byos, _parse_biopharma, _parse_dynamX, _parse_HDXWorkbench, \
        _parse_HDExaminer, _parse_peptide_list, _parse_custom, validate_and_parse_hxms_file

    demo.launch(debug=True,share=True)
