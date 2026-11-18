from dx_search.configs.common import PROJECT_SEARCH

CONFIG = {
    "log_level": "INFO",
    "output_file": "output/flagstat_files.json",
    "project_search": PROJECT_SEARCH,
    # Folder varies between runs (eggd_MultiQC/multiqc_data, <run>-multiqc_data, root), so don't filter on it
    "file_search": {
        "classname": "file",
        "name": r"^multiqc_samtools_flagstat\.txt$",
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
    # Removes ad-hoc copies such as /test_themis_multiQC/...
    "file_filters": [
        {
            "key": ["describe", "folder"],
            "mode": "regexp",
            "pattern": r"(?i)test",
            "outcome": "exclude",
        },
    ],
    "deduplicate": {
        "group_by": {
            "key": ["project"],
        },
    },
}
