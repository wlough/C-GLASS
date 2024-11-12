#!/bin/bash

# image_dir=$1
# if image_dir is not provided, use the current directory
if [ -z "$1" ]
then
    image_dir=$(pwd)
else
    image_dir=$1
fi

image_format="bmp"
image_prefix="test"
index_length=5
movie_name=image_prefix
movie_format="mp4"

frame_rate=200
frame_size="1080x720"
video_codec="libx264"
video_quality=25
pixel_format="yuv420p"


image_filename="${image_prefix}_%0${index_length}d.${image_format}"
movie_filename="${movie_name}.${movie_format}"


starting_dir=$(pwd)

run_command="ffmpeg"
ffmpegFLAGS=""
# overwrite output file without asking if it already exists
ffmpegFLAGS="$ffmpegFLAGS -y"
# frame rate (Hz)
ffmpegFLAGS="$ffmpegFLAGS -r $frame_rate"
# frame width x height (pixels)
ffmpegFLAGS="$ffmpegFLAGS -s $frame_size"
# input files path and format
ffmpegFLAGS="$ffmpegFLAGS -i $image_filename"
# video codec
ffmpegFLAGS="$ffmpegFLAGS -vcodec $video_codec"
# video quality, lower means better
ffmpegFLAGS="$ffmpegFLAGS -crf $video_quality"
# pixel format
ffmpegFLAGS="$ffmpegFLAGS -pix_fmt $pixel_format"
# output file
ffmpegFLAGS="$ffmpegFLAGS $movie_filename"


run_command="$run_command $ffmpegFLAGS"
echo $run_command
# Start the process
cd $image_dir
$run_command
cd $starting_dir
echo "Movie saved at ${image_dir}/$movie_path"
