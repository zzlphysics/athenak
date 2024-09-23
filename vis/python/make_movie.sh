#!/bin/bash

SIZE=50
DRICTION="z"
INPUT_FILE_PATH="/home/zzl/Project/athenak/output/mad/bin"
OUTPUT_FILE_PATH="/home/zzl/Project/athenak/output/mad"
OUTPUT_FILE_NAME="frames_log_dens_"$DRICTION"_"$SIZE

mpirun -n 8 python /home/zzl/Project/athenak/vis/python/make_movie.py -d $DRICTION -n log --r_max=$SIZE --vmin=1e-5 --vmax=1.5 --horizon_mask $INPUT_FILE_PATH dens $OUTPUT_FILE_PATH/$OUTPUT_FILE_NAME

# 设置帧率
FPS=30

# 检查 frames_log_rho 文件夹是否存在
if [ -d "$OUTPUT_FILE_PATH/$OUTPUT_FILE_NAME" ]; then
    echo "Encoding $OUTPUT_FILE_NAME.mp4"
    ffmpeg -hide_banner -loglevel error -y -r ${FPS} -f image2 -pattern_type glob -i "$OUTPUT_FILE_PATH/$OUTPUT_FILE_NAME/*.png" -vcodec libx264 -crf 22 -pix_fmt yuv420p "$OUTPUT_FILE_PATH/$OUTPUT_FILE_NAME.mp4"
else
    echo "frames_log_rho directory not found!"
fi
