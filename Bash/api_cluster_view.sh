#!/bin/bash
export PATH=$PATH:~/github/ProjectAPI/Bash
output_file_name=~/logs/viewoutput.csv

cd ~/outputs/

rm $output_file_name

ViewDateDirectory.sh /tigress/LEIFER/Mochi/APIData/ $1
ViewSortedStatus.sh