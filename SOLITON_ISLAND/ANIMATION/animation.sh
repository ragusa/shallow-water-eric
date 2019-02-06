ffmpeg -i 'movie.%03d.png' -vcodec h264 -filter:v "setpts=2.0*PTS" movie.mp4
rm *.vtk
rm *.png
