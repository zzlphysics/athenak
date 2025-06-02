#!/bin/bash

#SBATCH --job-name=plot
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=2
#SBATCH --time=1-00:00:00
#SBATCH --partition=debug
#SBATCH --output=plot_log/res_athk_plot_%j.out
#SBATCH --error=plot_log/res_athk_plot_%j.err

spack env activate athenak


SIZE=100
DRICTION="x"
VALUE="derived:beta_inv_rel"
# FILE_PATH="/home/zhangzelin/projects/athenak/output/20240929_hip_mad_a-05000_128-128-5_gamma14"
# FILE_PATH="/home/zhangzelin/projects/athenak/output/20240929_hip_mad_a-09375_128-128-5_gamma14"
# FILE_PATH="/home/zhangzelin/projects/athenak/output/20240929_hip_mad_a+00000_128-128-5_gamma14"
# FILE_PATH="/home/zhangzelin/projects/athenak/output/20240929_hip_mad_a+05000_128-128-5_gamma14"
# FILE_PATH="/home/zhangzelin/projects/athenak/output/20240929_hip_mad_a+09375_128-128-5_gamma14"
# FILE_PATH="/home/zhangzelin/projects/athenak/output/20240929_hip_mad_a+09800_128-128-5_gamma14"
# FILE_PATH="/home/zhangzelin/projects/athenak/output/20240929_hip_sane_a+05000_128-128-5_gamma13"
FILE_PATH="/home/zhangzelin/projects/athenak/output/20240929_hip_sane_a+09375_128-128-5_gamma13"

INPUT_FILE_PATH="$FILE_PATH/bin"
OUTPUT_FILE_PATH="$FILE_PATH"
OUTPUT_FILE_NAME="frames_log_"$VALUE"_"$DRICTION"_"$SIZE

srun -n 18 python /home/zhangzelin/projects/athenak/vis/python/make_movie.py -d $DRICTION -c jet -n log --r_max=$SIZE --vmin=1e-3 --vmax=1e3 --horizon_mask --notex $INPUT_FILE_PATH $VALUE $OUTPUT_FILE_PATH/$OUTPUT_FILE_NAME

# 设置帧率
FPS=30

# 检查文件夹是否存在
if [ -d "$OUTPUT_FILE_PATH/$OUTPUT_FILE_NAME" ]; then
    echo "Encoding $OUTPUT_FILE_NAME.mp4"
    ffmpeg -hide_banner -loglevel error -y -r ${FPS} -f image2 -pattern_type glob -i "$OUTPUT_FILE_PATH/$OUTPUT_FILE_NAME/*.png" -vcodec libx264 -crf 22 -pix_fmt yuv420p "$OUTPUT_FILE_PATH/$OUTPUT_FILE_NAME.mp4"
else
    echo "Directory not found!"
fi