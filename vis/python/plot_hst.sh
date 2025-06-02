#!/bin/bash

# 定义 filelist 和 valuelist
filelist=(
    "/home/zhangzelin/projects/athenak/output/20240929_hip_mad_a-05000_128-128-5_gamma14"
    "/home/zhangzelin/projects/athenak/output/20240929_hip_mad_a-09375_128-128-5_gamma14" 
    "/home/zhangzelin/projects/athenak/output/20240929_hip_mad_a+00000_128-128-5_gamma14" 
    "/home/zhangzelin/projects/athenak/output/20240929_hip_mad_a+05000_128-128-5_gamma14"
    "/home/zhangzelin/projects/athenak/output/20240929_hip_mad_a+09375_128-128-5_gamma14"
    "/home/zhangzelin/projects/athenak/output/20240929_hip_mad_a+09800_128-128-5_gamma14"
    "/home/zhangzelin/projects/athenak/output/20240929_hip_sane_a+05000_128-128-5_gamma13"
    "/home/zhangzelin/projects/athenak/output/20240929_hip_sane_a+09375_128-128-5_gamma13"
)
valuelist=("val1" "val2" "val3")

# 遍历 filelist
for file in "${filelist[@]}"
do
    # 遍历 valuelist
    for val in "${valuelist[@]}"
    do
        # 调用 python 脚本进行处理
        python plot_hst.py -i "$file/torus.user.hst" -o "$file/hst.png" -v "$val"
    done
done