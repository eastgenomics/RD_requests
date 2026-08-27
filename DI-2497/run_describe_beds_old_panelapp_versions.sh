#!/bin/bash
# Example run:$ bash run_describe_beds_old_panelapp_versions.sh R134 project-GQ5f2z04vFb3XzFvY181Yx63
TEST_CODE=$1
PROJECT_WITH_TEST_CODE=$2
PROJECT_NAME_WITH_TEST_CODE=$(dx describe --name "${PROJECT_WITH_TEST_CODE}")
panel_beds_R134=$(dx find data --name *"${TEST_CODE}"*.bed --project "${PROJECT_WITH_TEST_CODE}" --json | jq -r '.[] | .id') 

GENE_PANELS=$(dx cat project-Fkb6Gkj433GVVvj73J7x8KbV:file-J1jXFZj4XG7Qvj0PGZGg96Pg | awk '{ split($0, arr, "\t"); print arr[1]; }') 

# validate inputs 
# --- Check valid test code ---
if [[ $(grep -ic "${TEST_CODE}" <<< "${GENE_PANELS}") -eq 0 ]]
then
  echo "Error: '$TEST_CODE' is not a valid test code ."
  exit 1
fi

# --- Check valid  project for the test code  ---
# should be 002 and be a CEN/TWE project 

if [[ $(grep -ic  "002" <<< "${PROJECT_NAME_WITH_TEST_CODE}") -eq 0 ]]
then
  echo "Error: \"{$PROJECT_WITH_TEST_CODE}\"-\"${PROJECT_NAME_WITH_TEST_CODE}\" is not a valid 002 project ."
  exit 1
fi

if [[ $(grep -ic -E "CEN|TWE" <<< "${PROJECT_NAME_WITH_TEST_CODE}") -eq 0 ]]
then
  echo "Error: \"{$PROJECT_WITH_TEST_CODE}\"-\"${PROJECT_NAME_WITH_TEST_CODE}\" is not a valid 002 dias project ."
  exit 1
fi

for ID in ${panel_beds_R134[@]};do 
    echo "Dx describe: \"${ID}\""
    FILE_ID=$(dx describe "${ID}" --json | jq -r '.id') 
    FILE_NAME=$(dx describe "${ID}" --json | jq -r '.name') 
    PROJECT=$(dx describe "${ID}" --json | jq -r '.project')
    PROJECT_NAME=$(dx describe "${PROJECT}" --json | jq -r '.name')
    # get job id that generated the bed panel file 
    JOB_ID=$(dx describe "${ID}" --json | jq -r '.createdBy.job')
    # obtain job name, destination folder and analysis id from the job id 
    JOB_NAME=$(dx describe "${JOB_ID}" --json | jq -r '.name')
    OUT_FOLDER=$(dx describe "${JOB_ID}" --json | jq -r '.folder')
    ANALYSIS_ID=$(dx describe "${JOB_ID}" --json | jq -r '.analysis')
    # obtained anaysis workflow name 
    ANALYSIS_NAME=$(dx describe "${ANALYSIS_ID}" --json | jq -r '.name') 

    SAMPLE=$(cut -d '-' -f 2 <<< "${ANALYSIS_NAME}")
    echo "FILE_ID: \"${FILE_ID}\""
    echo "FILE_NAME: \"${FILE_NAME}\""
    echo "DESTINATION_FOLDER: \"${OUT_FOLDER}\""
    echo "PROJECT: \"${PROJECT_NAME}\" (\"${PROJECT}\")"
    echo "JOB_NAME: \"${JOB_NAME}\""
    echo "ANALYSIS: \"${ANALYSIS_NAME}\" (\"${ANALYSIS_ID}\")"
    echo "SAMPLE: \"${SAMPLE}\""
    echo "------------------End of processing ------------------"
    done