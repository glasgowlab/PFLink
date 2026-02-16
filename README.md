---
title: PFLink
emoji: 📊
colorFrom: green
colorTo: blue
sdk: gradio
sdk_version: 6.5.1
app_file: app.py
pinned: false
license: mit
---
Click here to use the [website](https://huggingface.co/spaces/glasgow-lab/PFLink).

## Files needed to run PFLink

* **HDX Workbench**
    * The default file export.
* **HDExaminer**
    * Two files are needed: "pool peptide results" (CSV) and "pool spectra export" (folder).
* **DynamX**
    * The default file export.
* **BioPharma Finder**
    * HDXData.csv
* **Custom**
    * The default custom file is found in the `test_data` directory of this repository.

---

## Test data files

Example datasets are available in the [**`test_data`**](https://huggingface.co/spaces/glasgow-lab/PFLink/resolve/main/test_data/test_data.zip?download=true) directory in this repository.

All ecDHFR data are from [Pigeon-Feather](https://doi.org/10.1101/2024.08.04.606547).

* `SpecExport.zip` contains the isotopic mass envelope data for the ecDHFR dataset. This was generated from HDExaminer (pool spectra export).
* `ecDHFR_tutorial.csv` contains the isotopic uptake data for the ecDHFR dataset. This was also generated from HDExaminer (pool peptide results).
* `flag-DHFR.txt` contains all the protein and experimental metadata needed to generate the HXMS file.
* `rangeslist.csv` contains the peptides that should be kept when using PFLink.
* `custom_empty.csv` is an empty PFLink file that a user can fill in if they are not using a supported program.
* `custom.csv` is a filled example of the custom format with fake data values, as an example.

---

## Installing and running PFLink locally

### Installing via Conda
```
git clone https://huggingface.co/spaces/glasgow-lab/PFLink
cd PFLink
conda create -n PFLink python=3.11
conda activate PFLink
pip install gradio
pip install git+https://github.com/glasgowlab/PIGEON-FEATHER.git
pip install numba
pip install pyopenms
pip install hdxrate
pip install rich-argparse
```

### Running the website locally via Conda
```
conda activate PFLink
cd PFLink # if not already in the directory
python3 app.py
```

### Running the website locally via Docker
If you prefer running PFLink from the Docker rather than the local Conda installation, skip the instructions above and do this instead:

```
docker run -it -p 7860:7860 --platform=linux/amd64

registry.hf.space/glasgow-lab-pflink:latest python app.py
```

### Running the scripts locally via Conda
```
conda activate PFLink
cd PFLink # if not already in the directory
python3 pflink_isotopic_uptake.py # for non-envelope data processing
python3 pflink_envelope.py # for envelope data processing
python3 pflink_merge.py # for merging different .HXMS files
```

---

## Using the isotopic_uptake data processing script

This script is meant to convert DynamX, HDX Workbench, HDExaminer, BioPharma, or our custom export format into an `.HXMS` file. This is intended for centroid data only.

**Arguments:**
```
--help
Show this help message and exit.

--input_csv_path INPUT_CSV_PATH
Path to the file containing the raw input CSV/data file.

--output_hxms_path OUTPUT_HXMS_PATH
Path to the output directory where the final .HXMS files will be saved.

--flags_file_path FLAGS_FILE_PATH
Path to a flags or configuration file to override default settings.

--peptide_list PEPTIDE_LIST
Path to a CSV/TXT file containing a list of peptides to filter or process.

--saturation SATURATION
D2O saturation percentage (0.0 to 1.0) to record in the file header.

--ph PH
The measured pH (pH Read) value of the exchange buffer.

--temperature TEMPERATURE
The temperature of the exchange reaction in Kelvin (K).

--protein_name PROTEIN_NAME
The name of the protein.

--protein_state PROTEIN_STATE
The experimental state of the protein (e.g., 'Apo', 'Bound', 'Mutant').

--protein_sequence PROTEIN_SEQUENCE
The full amino acid sequence of the protein.

--include_exclude {include,exclude}
Determines if the peptide list file should be used to include or exclude peptides from processing.

--file_type {DynamX,HDXworkbench,HDExaminer,BioPharma,Custom}
Input CSV format type

**`input_csv_path`** and **`output_hxms_path`** are strictly required. However, you must either provide a complete flags file or use all options. The only two optional parameters are `include_exclude` and `peptide_list`.
```
**Command usage:**
```
python3 pflink_isotopic_uptake.py --input_csv_path {path} --output_hxms_path {path} --flags_file_path {path}
```
or
```
python3 pflink_isotopic_uptake.py --input_csv_path {path} --output_hxms_path {path} --saturation {saturation} --ph {ph} --protein_name {name} --protein_state {state} --protein_sequence {sequence} --file_type {file_type}
```
Examples of the flags file and peptide list can be found in the [`test_data`](https://huggingface.co/spaces/glasgow-lab/PFLink/resolve/main/test_data/test_data.zip?download=true) folder.

---

## Using the isotopic envelope data processing script

This script is meant to convert the HDExaminer export format into an `.HXMS` file. This is intended for HDExaminer envelope data only.

**Arguments:**
```
--help
Show this help message and exit.

--input_csv_path INPUT_CSV_PATH
Path to the file containing the raw input CSV/data file.

--output_hxms_path OUTPUT_HXMS_PATH
Path to the output directory where the final .HXMS files will be saved.

--flags_file_path FLAGS_FILE_PATH
Path to a flags or configuration file to override default settings.

--peptide_list PEPTIDE_LIST
Path to a CSV/TXT file containing a list of peptides to filter or process.

--raw_spectra_path RAW_SPECTRA_PATH
Path of the HDExaminer spectra folder.

--saturation SATURATION
D2O saturation percentage (0.0 to 1.0) to record in the file header.

--ph PH
The measured pH (pH Read) value of the exchange buffer.

--temperature TEMPERATURE
The temperature of the exchange reaction in Kelvin (K).

--protein_name PROTEIN_NAME
The name of the protein.

--protein_state PROTEIN_STATE
The experimental state of the protein (e.g., 'Apo', 'Bound', 'Mutant').

--protein_sequence PROTEIN_SEQUENCE
The full amino acid sequence of the protein.

--include_exclude {include,exclude}
Determines if the peptide list file should be used to include or exclude peptides from processing.

--save_match
Whether to save the match data. Default is False.

--save_fine_match
Whether to save the fine match data (uncentroided, fine structure). Default is False.

--file_type {HDExaminer}
Input CSV format type
```

**`input_csv_path`**, **`output_hxms_path`**, and **`raw_spectra_path`** are strictly required. However, you must either provide a complete flags file or use all options. The only two optional parameters are `include_exclude` and `peptide_list`.

**Command usage:**
```
python3 pflink_envelope.py --input_csv_path {path} --output_hxms_path {path} --raw_spectra_path {path} --flags_file_path {path}
```
or
```
python3 pflink_isotopic_uptake.py --input_csv_path {path} --output_hxms_path {path} --raw_spectra_path {path} --saturation {saturation} --ph {ph} --protein_name {name} --protein_state {state} --protein_sequence {sequence} --file_type {file_type}
```
Examples of the flag file and peptide list can be found in the [`test_data`](https://huggingface.co/spaces/glasgow-lab/PFLink/resolve/main/test_data/test_data.zip?download=true) folder.

---

## Using the HXMS merge script

This script is for combining replicates -- different experimental `.HXMS` runs of the same protein or state.

**Arguments:**
```
--help
Show this help message and exit.

--input_hxms_path1 INPUT_HXMS_PATH1
Path to the first file containing the .HXMS data.

--input_hxms_path2 INPUT_HXMS_PATH2
Path to the second file containing the .HXMS data.

--output_hxms_path OUTPUT_HXMS_PATH
Path to the output directory where the final .HXMS files will be saved.
```

All options are required to run the script.

**Command usage:**
```
python3 pflink_merge.py --input_hxms_path1 {path} --input_hxms_path2 {path} --output_hxms_path {path}
```

---

## Using the website version on HuggingFace

### Example datasets

* Use the **"Isotopic uptake data load example"** to load the ecDHFR dataset for isotopic uptake only.
* Use the **"Envelope data load example"** to load the ecDHFR dataset for both isotopic uptake and envelopes.

Once selected, all fields will be autofilled with the correct experimental and protein information.

Click the **"Process Data"** button at the bottom to generate the HXMS output file.