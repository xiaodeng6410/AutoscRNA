#!/bin/bash

CALL_DIR=$1
cat ${CALL_DIR}/filename.txt | while read id; do
    /home/dengys/anaconda3/envs/cellbender/bin/python /data_result/dengys/scRNA/AutoScRNA/Script/CellbendertoMartix.py ${CALL_DIR}/cellbender_out/${id}
done

