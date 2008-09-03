#!/bin/bash
# 
# This will make a movie from PNG's in a folder and name the movie after the foler name
# 

for i in *; do
	if [ -d "${i}"]; then
		echo -n "Processing '${i}'\t\t"
    mencoder "mf://${i}/*.png" -mf fps=5 -o ${i}.avi -ovc lavc -lavcopts vcodec=mpeg4 &> /dev/null && echo "[DONE]" || (echo "[ !! ]" && exit );
  fi
done;
