def raise_if_archived(files, description):
    """
    Raise if any DNAnexus file is not live, listing every such file.

    Parameters
    ----------
    files : list
        File dicts from dx_search.search()["files"], searched with archivalState described
    description : str
        What the files are, for the error message
    """
    states = [(f["project"], f["id"], f["describe"].get("archivalState")) for f in files]
    not_live = [item for item in states if item[2] != "live"]
    if not_live:
        listing = "\n".join(f"  {project} {file_id} ({state})" for project, file_id, state in not_live)
        raise RuntimeError(f"{len(not_live)} {description} files are not live and need unarchiving:\n{listing}")
