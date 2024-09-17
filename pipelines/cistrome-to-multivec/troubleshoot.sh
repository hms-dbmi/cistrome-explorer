#!/bin/bash

# This script can be used to detect badly-formatted JSON (typically HTML server error pages)
# downloads of the .data.json and .metadata.json files,
# and delete any bad files.

# bash ./troubleshoot.sh .
# bash ./troubleshoot.sh /n/scratch3/users/m/mk596/cistrome-explorer

BASE_DIR=$1

echo "Running in $BASE_DIR"

for f in $BASE_DIR/data/raw/*.data.json; do
    if ! cat $f | jq '.[0].url' -r; then
        rm $f
    fi
done