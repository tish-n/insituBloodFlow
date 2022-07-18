#!/bin/sh
# this is a transfer script for output in WSL
# replace the directory below with whichever directory you want to move the output to
dynamicBoundary=$PWD

echo 

cd /mnt/d/PViewData

rm -rf *

cp $dynamicBoundary/tmp/* /mnt/d/PViewData

cd ..
