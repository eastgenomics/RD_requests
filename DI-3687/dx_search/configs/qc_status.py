from dx_search.configs.common import PROJECT_SEARCH

CONFIG = {
    "log_level": "INFO",
    "output_file": "output/qc_status_files.json",
    "project_search": PROJECT_SEARCH,
    # <run>_QC_Status.xlsx sits in the root of each b38 project
    "file_search": {
        "classname": "file",
        "name": r"(?i)_QC_Status\.xlsx$",
        "name_mode": "regexp",
        "describe": {
            "folder": True,
            "name": True,
            "created": True,
            "modified": True,
            "size": True,
            "archivalState": True,
        },
    },
    "deduplicate": {
        "group_by": {
            "key": ["project"],
        },
    },
}
