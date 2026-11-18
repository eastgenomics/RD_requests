import fnmatch
import json
import logging
import operator
import re
from collections import defaultdict

logger = logging.getLogger(__name__)


def _traverse_dict(d, path):
    for key in path:
        try:
            d = d[key]
        except KeyError as e:
            raise KeyError(
                f"Error when traversing dictionary: key '{key}' not found in the provided dictionary {json.dumps(d)}."
            ) from e
        # If a value along the specified path is None, raise rather than treat it as a non-match
        if d is None:
            raise ValueError(f"No value found for key '{key}' in path {path}.")
    return d


def _numeric_match(value, expression):
    ops = {
        ">": operator.gt,
        ">=": operator.ge,
        "<": operator.lt,
        "<=": operator.le,
        "==": operator.eq,
        "!=": operator.ne,
    }

    expression = str(expression).strip()

    for op in sorted(ops, key=len, reverse=True):
        if expression.startswith(op):
            target = float(expression[len(op) :])
            return ops[op](float(value), target)

    raise ValueError(f"Invalid numeric expression: {expression}")


def _matches(dx_dict, key, mode, pattern, outcome="include"):
    if mode not in ("glob", "bool", "regexp", "exact", "numeric", "list"):
        raise ValueError(f"Unknown filter mode: {mode}")

    value = _traverse_dict(dx_dict, key)

    if mode == "glob":
        if "folder" in key:
            value = value.strip("/").lower()
            pattern = pattern.strip("/").lower()
            matched = fnmatch.fnmatch(value, pattern)
        else:
            matched = fnmatch.fnmatchcase(str(value), pattern)

    if mode == "bool":
        matched = bool(value) == bool(pattern)

    if mode == "regexp":
        matched = re.search(pattern, value) is not None

    if mode == "exact":
        matched = value == pattern

    if mode == "numeric":
        matched = _numeric_match(value, pattern)

    if mode == "list":
        matched = value in pattern

    return not matched if outcome == "exclude" else matched


def filter_dx_dicts(dx_dicts, filters):
    results_dict = {
        "included": dx_dicts.copy(),
        "excluded": {},
    }
    for dx_id, dx_dict in dx_dicts.items():
        for filter_config in filters:
            matched = _matches(dx_dict=dx_dict, **filter_config)

            if not matched:
                del results_dict["included"][dx_id]
                results_dict["excluded"][dx_id] = dx_dict
                break

    return results_dict


def _get_newest(file_group):
    newest_idx = None
    for idx in range(len(file_group)):
        if newest_idx is None:
            newest_idx = idx
            continue

        cur_create_ts = file_group[idx]["describe"]["created"]
        newest_create_ts = file_group[newest_idx]["describe"]["created"]

        if cur_create_ts > newest_create_ts:
            newest_idx = idx
        elif cur_create_ts == newest_create_ts:
            newest_mod_ts = file_group[newest_idx]["describe"]["modified"]
            cur_mod_ts = file_group[idx]["describe"]["modified"]
            if cur_mod_ts > newest_mod_ts:
                newest_idx = idx
            if cur_mod_ts == newest_mod_ts:
                logger.warning(
                    "Files %s and %s have the same created and modified timestamps. "
                    "Arbitrarily selecting one as the newest.",
                    file_group[newest_idx]["describe"]["name"],
                    file_group[idx]["describe"]["name"],
                )
    newest = file_group.pop(newest_idx)
    not_newest = file_group
    return newest, not_newest


def deduplicate_files(files, group_by):
    file_groups = defaultdict(list)
    for file in files:
        group_by_key = _traverse_dict(file, group_by["key"])
        if group_by.get("substring"):
            if not isinstance(group_by_key, str):
                raise ValueError(f"Group key must be a string to apply substring pattern, got {type(group_by_key)}")
            substring_pattern = group_by["substring"]["pattern"]
            substring_match = re.search(substring_pattern, group_by_key)
            if not substring_match:
                raise ValueError(
                    f"Specified substring pattern '{substring_pattern}' could not be found in group by key "
                    f"'{group_by_key}'"
                )
            group_by_key = substring_match.group(1)
        file["group_by_key"] = group_by_key
        file_groups[group_by_key].append(file)

    result_dict = {"deduplicated_files": [], "duplicate_files": []}
    for file_group in file_groups.values():
        if len(file_group) == 1:
            result_dict["deduplicated_files"].append(file_group[0])
            continue
        newest, not_newest = _get_newest(file_group)
        result_dict["deduplicated_files"].append(newest)
        result_dict["duplicate_files"].extend(not_newest)
    return result_dict
