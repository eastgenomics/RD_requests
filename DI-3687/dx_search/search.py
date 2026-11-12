import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Any

from .config import configure_logging, load_config
from .dx_filter import deduplicate_files, filter_dx_dicts
from .find import find_files_in_project, find_projects

logger = logging.getLogger(__name__)


def search(config: dict[str, Any]) -> dict[str, list]:
    """
    Find projects and files in DNAnexus according to a config dict.

    Parameters
    ----------
    config : dict
        Search config (see configs/flagstat.py for the expected keys)

    Returns
    -------
    dict
        Keys: files, excluded_files, duplicates, excluded_projects. Each file
        is the dxpy find_data_objects result dict (id, project, describe).
    """
    output_dict = {
        "files": [],
        "excluded_files": [],
        "duplicates": [],
        "excluded_projects": [],
    }

    logger.info("Finding projects...")
    projects = find_projects(config["project_search"])

    if not projects:
        logger.info(
            "No projects found matching the search criteria. Search criteria:\n%s",
            json.dumps(config["project_search"], indent=2),
        )
        return output_dict

    project_filter_configs = config.get("project_filters")
    if project_filter_configs:
        logger.info("Filtering projects...")
        project_filter_results = filter_dx_dicts(projects, project_filter_configs)
        output_dict["excluded_projects"].extend(project_filter_results["excluded"].values())
        projects = project_filter_results["included"]

        if not projects:
            logger.info(
                "No projects found after filtering. \nProject filters:\n%s",
                json.dumps(project_filter_configs, indent=2),
            )
            return output_dict
    logger.info("Found %d projects matching criteria.", len(projects))

    logger.info("Finding files in projects...")
    file_filter_configs = config.get("file_filters")
    for project_id in projects:
        files = find_files_in_project(project_id, config["file_search"])

        if file_filter_configs:
            file_filter_results = filter_dx_dicts(files, file_filter_configs)
            output_dict["excluded_files"].extend(file_filter_results["excluded"].values())
            output_dict["files"].extend(file_filter_results["included"].values())
        else:
            output_dict["files"].extend(files.values())

    if not output_dict["files"]:
        logger.info("No files found.")
        return output_dict

    deduplicate_config = config.get("deduplicate")
    if deduplicate_config:
        logger.info("Deduplicating files globally...")
        result = deduplicate_files(files=output_dict["files"], group_by=deduplicate_config["group_by"])
        output_dict["files"] = result["deduplicated_files"]
        output_dict["duplicates"] = result["duplicate_files"]

    return output_dict


def parse_arguments(argv=None):
    parser = argparse.ArgumentParser(description="Search for files in DNAnexus.")
    parser.add_argument("-c", "--config-file", help="Path to config file", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_arguments(argv)
    config = load_config(args.config_file)
    configure_logging(config.get("log_level", "INFO"))
    logger.info("Loaded configuration from\n %s", json.dumps(config, indent=2))

    output_dict = search(config)

    output_file = config.get("output_file")
    if output_file:
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        with open(output_file, "w") as f:
            json.dump(output_dict, f, indent=2)
        logger.info("Output saved to %s", output_file)

    if config.get("emit_json") is True:
        sys.stdout.write(json.dumps(output_dict, indent=2))


if __name__ == "__main__":
    main()
