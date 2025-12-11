#!/bin/bash

SIZE=50
DRICTION="z"

# VALUE="dens"
# VMIN=1e-5
# VMAX=15
# CMAP="viridis"

# VALUE="derived:beta_inv_rel"
# VMIN=1e-3
# VMAX=1e3
# CMAP="jet"
# CMAP="viridis"

VALUE="derived:sigma_rel"
VMIN=1e-8
VMAX=1e6
# CMAP="jet"
CMAP="viridis"


FILE_PATH="/home/zzl/kx-4t/athenak/output/20240929_cuda_mad_a09375_128-128-128-8"
INPUT_FILE_PATH="$FILE_PATH/bin"
OUTPUT_FILE_PATH="$FILE_PATH"
OUTPUT_FILE_NAME="frames_log_"$VALUE"_"$DRICTION"_"$SIZE

mpirun -n 6 python ./make_movie.py -d $DRICTION -n log -c $CMAP --r_max=$SIZE --vmin=$VMIN --vmax=$VMAX --horizon_mask --notex $INPUT_FILE_PATH $VALUE $OUTPUT_FILE_PATH/$OUTPUT_FILE_NAME

# 设置帧率
FPS=30

# 检查 frames_log_rho 文件夹是否存在
if [ -d "$OUTPUT_FILE_PATH/$OUTPUT_FILE_NAME" ]; then
    echo "Encoding $OUTPUT_FILE_NAME.mp4"
    ffmpeg -hide_banner -loglevel error -y -r ${FPS} -f image2 -pattern_type glob -i "$OUTPUT_FILE_PATH/$OUTPUT_FILE_NAME/*.png" -vcodec libx264 -crf 22 -pix_fmt yuv420p "$OUTPUT_FILE_PATH/$OUTPUT_FILE_NAME.mp4"
else
    echo "frames_log_rho directory not found!"
fi
