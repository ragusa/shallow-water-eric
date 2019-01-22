rm *.plt
cp ../bath.plt .
mv ../hpz_movie_* .
gnuplot -p movie.gnu
ffmpeg -i 'movie.%03d.png' -vcodec h264 -filter:v "setpts=2.0*PTS" movie.mp4
#ffmpeg -i $1.%4d.png -codec:v libtheora -qscale:v 7 -filter:v "setpts=2.5*PTS" $1-vid.ogv
rm *.png