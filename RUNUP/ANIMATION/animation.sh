#ffmpeg -i $1.%4d.png -vcodec h264 -filter:v "setpts=2.5*PTS" $1-vid.ogv
ffmpeg -i $1.%4d.png -codec:v libtheora -qscale:v 8 -filter:v "setpts=2.5*PTS" $1-vid.ogv
