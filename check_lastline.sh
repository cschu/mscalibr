#!/bin/bash




for fn in $(find . -type f)
do
  if tail -n 1 $fn | grep --quiet "ION ENDS" $fn; then
    echo $fn "ION ENDS PRESENT"
  else
    echo $fn
  fi
done
