import argparse
import json
import re

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

URL = "https://panelapp.genomicsengland.co.uk/api/v1/panels/signedoff/"


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Write the ID, name and R code(s) of every signed-off PanelApp "
            "panel to a JSON file"
        )
    )
    parser.add_argument(
        "output", help="Path of the JSON file to write, e.g. panels.json"
    )
    return parser.parse_args()


def get_signedoff_panels():
    """Return {panel_id: {"rcode": ..., "panel_name": ...}} for all
    signed-off panels"""
    with requests.Session() as session:
        session.headers["Accept"] = "application/json"
        session.mount(
            "https://",
            HTTPAdapter(
                max_retries=Retry(
                    total=6,
                    backoff_factor=1,
                    status_forcelist=[429, 500, 502, 503, 504],
                )
            ),
        )

        # Get name and relevant disorders
        p_dict = {}
        url = URL
        while url:
            resp = session.get(url, timeout=30)
            resp.raise_for_status()
            data = resp.json()
            for panel in data["results"]:
                rcodes = [
                    x
                    for x in panel["relevant_disorders"]
                    if re.match(r"R\d+", x)
                ]
                p_dict[panel["id"]] = {
                    "rcode": ",".join(rcodes),
                    "panel_name": panel["name"],
                }
            url = data["next"]

    return p_dict


def main():
    args = parse_args()
    p_dict = get_signedoff_panels()

    with open(args.output, "w") as f:
        json.dump(p_dict, f)

    print(f"Wrote {len(p_dict)} panels to {args.output}")


if __name__ == "__main__":
    main()