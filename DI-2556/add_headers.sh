#!/bin/bash
# save comparison output
touch difference_between_old_and_new_static_bed_panel.txt
# track number of static beds processed
count=0

# define old static bed file path and version number
old_static_bed_path="old_static_beds/v1.0.0*.bed"
old_version=v1.0.0

# define new static bed file path and version number 
new_static_bed_path="new_static_beds"
new_version=v1.0.1


#update old static bed file by adding a header detailing the updated versions and save in new static file folder.
# Then use diff to compare old and new static bed files 
for old_panel_bed in $old_static_bed_path ; do
  let count++
  printf  "Number of static beds being processed: %d\n" $count

  # get filename and update the version to v1.0.1 in the filename and make file path 
  new_panel_bed=$(basename $old_panel_bed | sed "s/$old_version/$new_version/") 
  new_panel_bed_path="${new_static_bed_path}/${new_panel_bed}"
  printf  "File: %s to be updated to ----> %s ......\n" $old_panel_bed $new_panel_bed

  # add header and save new static bed file
  header="#Panel_BED_version:$new_version"
  sed  "1i $header" $old_panel_bed >  $new_panel_bed_path
  printf  "New bed panel file path,  %s saved \n"  $new_panel_bed_path
  
  #Save diff result in a file 
  diff_output=$(diff -s $old_panel_bed $new_panel_bed_path)
  printf "Difference between $old_panel_bed and $new_panel_bed_path: \n$diff_output\n" >> difference_between_old_and_new_static_bed_panel.txt

done