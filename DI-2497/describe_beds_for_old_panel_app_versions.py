import requests
from pprint import pprint
import argparse
import dxpy
import pandas as pd 
from datetime import datetime
# -name *R134*.bed --created-after="2020-02-17" --all-projects --created-before="2023-03-22"
#s$ python describe_beds_for_old_panel_app_versions.py --test_code "R134.1" --created_after "2020-02-17"  --created_before "2023-03-22"

def describe_projects(list_project_id):
    count=0
    list_projects_with_test_code=[]
    for proj_id in list_project_id:
        count+=1
        name=dxpy.describe(proj_id)["name"]
        proj_id=dxpy.describe(proj_id)["id"]
        proj_type=name.split("_")[0]
        date_created=name.split("_")[1]
        list_projects_with_test_code.append({   
            "project_number":count,
            "name":name,
            "proj_id": proj_id, 
            "project_type":proj_type,
            "date_created":date_created
        })

    
    return list_projects_with_test_code


                 

    


def find_production_projects(test_code, start, end):
    """
    Find files in DNAnexus project by name

    Parameters
    ----------
    project_name : str
        project name
    start : str (required)
        start date to look for projects from
    end: str (requered)
        end date to look for projects until

    Returns
    -------
    projects : list
        list of DNAnexus projects
    """
    analyses=test_code.split(".")[1]
    test_code=test_code.split(".")[0]
    
    print(test_code)
    file_with_test_code="".join(["*",test_code,"*",".bed"])
    print(file_with_test_code)
    files = list(
        dxpy.find_data_objects(
            name=file_with_test_code,
            name_mode="glob",
            created_after=start,
            created_before=end,
            describe=True

        )
    )
    #obtain string ID of the project in which the result was found.
    print(f"Number of panel beds with test code {test_code}: {len(files)}")
    project_ids=[file['project']for file in files]
    
    print(f"Number of projects with test code {test_code}: {len(project_ids)}")
    print(f"Number of unique projects with test code {test_code}: {len(list(set(project_ids)))}")
    project_list=describe_projects(list(set(project_ids)))

    prod_projects=[proj for proj in project_list if proj["project_type"] == "002"]
    print(f"Number of unique production projects with test code {test_code}: {len(prod_projects)}")
    return prod_projects



def main():
    """Entry point."""
    parser = argparse.ArgumentParser(
        description="Get signed off panel for historical versions to locate old panel beds"
    )
    parser.add_argument(
        '--test_code', required=True,
        help="String for timestamp e.g. 2020-02-17"
    )
    parser.add_argument(
        '--created_after', required=True,
        help=""
    )

    parser.add_argument(
        '--created_before', required=True,
        help=""
    )
    args = parser.parse_args()

    print(f"Obtain list of projects which has test_code : {args.test_code}")
    find_production_projects(test_code=args.test_code, start=args.created_after, end=args.created_before)

    #print(f"Number of projects with test code {args.test_code} created after({args.created_after}) and before ({args.created_before})  : len({lists_projects_with_test_code})")


    
if __name__ == '__main__':
    main()

