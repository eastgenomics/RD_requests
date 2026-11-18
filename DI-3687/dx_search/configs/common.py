# Shared by flagstat.py and qc_status.py so both searches cover the same projects
PROJECT_SEARCH = {
    "name": r"^002_26.*_38_CEN$",
    "name_mode": "regexp",
    "created_after": "2026-07-01",
    "created_before": "2026-10-30",
    "describe": {
        "name": True,
    },
}
