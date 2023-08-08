#!/bin/bash
export PATH=$PATH:~/github/ProjectAPI/Bash

#start_point=-1 #redo everything
start_point=$(ScriptToOrdering.sh max) #continue analysis
#start_point=$(($start_point + 1))
#start_point=3

echo $start_point

cd ~/outputs/

InitalizeLogs.sh

ProcessDateDirectory.sh /tigress/LEIFER/Mochi/APIData/ $start_point