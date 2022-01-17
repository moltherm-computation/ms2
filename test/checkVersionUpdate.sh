#!/bin/bash

# Get difference of HEAD and previous commits
# The sqash commit does not count, therefore

diffOutput1=$(git diff HEAD~1 HEAD -- src/ms2_version.F90)

if [ -z "$diffOutput1" ]; # empty string = no update
then
    echo "Version not updated."
    exit 1
fi
