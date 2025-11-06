#!/bin/bash
# save comparison output
diff_output_file="$(date +'%Y%m%d')"_difference_between_old_and_new_static_bed_panel.txt
touch "${diff_output_file}"
# track number of static beds processed
count=0

# define old static bed file path and version number
old_static_bed_path=$1
old_version=$2


# define new static bed file path and version number 
new_static_bed_path="/home/dnanexus/new_static_beds"
new_version=$3
mkdir -p "${new_static_bed_path}" 

#update old static bed file by adding a header detailing the updated versions and save in new static file folder.
# Then use diff to compare old and new static bed files 
for old_panel_bed_path in ${old_static_bed_path}* ; do
  let count++
  printf  "Number of static beds being processed: %d\n" "$count"

  # get filename and update the version to v1.0.1 in the filename and make file path for new static bed 
  old_basename=$(basename "$old_panel_bed_path")
  new_panel_bed="${old_basename/$old_version/$new_version}"
  new_panel_bed_path="${new_static_bed_path}/${new_panel_bed}"
  printf  "File: %s to be updated to ----> %s ......\n" "$old_basename" "$new_panel_bed"

  # add header and save new static bed file
  header="#assembly=GRCh38,version=${new_version}"
  if sed  "1i ${header}" "${old_panel_bed_path}" >  "${new_panel_bed_path}"; then 
    printf  "New bed panel file path,  %s saved \n"  "$new_panel_bed_path"
    #Save diff result in a file 
    diff_output=$(diff -s "${old_panel_bed_path}" "${new_panel_bed_path}")
    printf "Difference between %s and %s :\n%s\n" "$old_panel_bed_path" "$new_panel_bed_path" "$diff_output" >> "$diff_output_file"
  else
    printf "Error: Failed to process %s\n" "$old_panel_bed_path" >&2
  fi 

done