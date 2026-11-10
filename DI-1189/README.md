## Create PanelApp JSON file

### Description
This script is used to create a PanelApp JSON file which is used as an input to the [eggd_optimised_filtering](https://github.com/eastgenomics/eggd_optimised_filtering) app to annotate the MOI of a gene in the context of a panel.

### Setup
Create a virtual environment and install the requirements:

```bash
python3 -m venv ~/venvs/panel_updates
source ~/venvs/panel_updates/bin/activate
pip install -r requirements.txt
```

Tested with Python 3.14. The pinned `pandas` 3.x requires Python 3.11 or later.

Note: the `panelapp` package does not declare `requests` as a dependency, so it is listed explicitly in `requirements.txt`.

### Inputs
The script `create_panelapp_dump.py` takes inputs:
- `--genepanels (str)`: genepanels file (optional)
- `--extra_panels (str)`: Any extra panel IDs to retain, comma separated (optional, requires `--genepanels` to be given)
- `--outfile_name (str)`: Name of output PanelApp JSON file

### Running
Example command:
`python create_panelapp_dump.py --genepanels 241024_genepanels.tsv --extra_panels 1570 --outfile_name 241030_panelapp_dump.json`

### How it works
The script:
- Gets the list of signed-off panels from PanelApp using the panelapp library.
- If provided, keeps only panels with an ID present in the input `--genepanels` file, plus any IDs given as `--extra_panels`, and fetches full data only for those panels. A warning is printed for any requested panel ID that is not signed off in PanelApp.
- Formats these panels as dictionaries and adds to list
- Checks for any duplicate genes or regions (by gene symbol/region name) within a panel
  - If only the MOI is different between any duplicates, retains the first instance with the MOI replaced to 'Other'.

#### PanelApp rate limiting
The PanelApp API rate-limits clients (HTTP 429), and the panelapp library gives up on the first error response. To handle this, the script waits at least 0.5 seconds between API calls and retries on HTTP 429 and 5xx responses up to 8 times with exponential backoff, following any `Retry-After` header from PanelApp. While it is being rate-limited, the script may pause for a few minutes on a single panel; this is expected and not a hang.

If any panel still cannot be fetched after retrying, the script exits with an error listing the failed panels and **no output file is written**. If this happens, wait a few minutes and run it again.

### Output
The output of the script is a JSON file named by `--outfile_name`, containing an array of panels.
