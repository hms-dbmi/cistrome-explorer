#!/usr/bin/env bash

cd "$(dirname "$0")"

HIGLASS_SERVER_PATH=$1

echo "Using higlass-server located at ${HIGLASS_SERVER_PATH}"

for D in $(ls -1 ./data/processed | sort); do
    if [ ${D: -4} == ".mv5" ]; then
        echo "Ingesting ${D}"
        python $HIGLASS_SERVER_PATH/manage.py ingest_tileset \
            --uid $D \
            --filename ./data/processed/$D \
            --filetype multivec \
            --datatype matrix
    fi
done