#!/bin/bash
CALL_DIR=$1 
cat ${CALL_DIR}/filename.txt | while read id;do 
  CUDA_VISIBLE_DEVICES=1 /home/dengys/anaconda3/envs/scvi2/bin/python /data_result/dengys/scRNA/AutoScRNA/Script/scanvi_annotate_10x.py ${CALL_DIR}/cellbender_out/${id}/10x_mtx 
done



