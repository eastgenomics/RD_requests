import importlib.util
import logging
import sys
from pathlib import Path
from typing import Any


def configure_logging(log_level: str) -> None:
    handlers = []
    log_format = "%(asctime)s - %(levelname)s - %(name)s - %(message)s"

    file_handler = logging.FileHandler("dx_search.log", mode="a")
    file_handler.setFormatter(logging.Formatter(log_format))
    handlers.append(file_handler)

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logging.Formatter(log_format))
    handlers.append(console_handler)

    logging.basicConfig(level=log_level, handlers=handlers)


def load_config(config_path: Path) -> dict[str, Any]:
    """
    Load the configuration file

    Parameters
    ----------
    config_path : Path
        Path to the Python configuration file

    Returns
    -------
    dict
        Configuration dictionary loaded from the specified file
    """
    spec = importlib.util.spec_from_file_location(config_path.name.replace(".py", ""), config_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[config_path.name] = module
    spec.loader.exec_module(module)

    if not hasattr(module, "CONFIG"):
        raise ValueError(f"Config file {config_path} must define a CONFIG object")

    return module.CONFIG
