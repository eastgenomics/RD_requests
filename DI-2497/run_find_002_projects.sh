#!/bin/bash

dx find data --name *R134*.bed --created-after="2020-02-17" --all-projects --created-before="2023-03-22" --json | jq -r '.[].project'| uniq > projects_with_R134.1_old_version_v1.2.txt

touch projects_002_with_R134.1_old_version_v1.2.txt
while read project; do
    PROJ_NAME=$(dx describe "$project" --json | jq '.name')
    if grep -q "002" <<< $PROJ_NAME 
    then 
        dx describe "$project" --json | jq -r '[.name, .id] | @csv'  >> projects_002_with_R134.1_old_version_v1.2.txt
    else
        echo "Not a 002 project: $PROJ_NAME" 
    fi
done < projects_with_R134.1_old_version_v1.2.txt
