#!/bin/sh

#  rectangle.sh
#  Created by Eric Tovar on 10/10/2018.
#
source ~/.bash_profile

# input domain D = [x1,x2] x [y1,y2] and refinment level 3,4,5,etc
# as x1 x2 y1 y2 refinment
x1="$1"
x2="$2"
y1="$3"
y2="$4"

# edit poly file
# sed -i -e 's/abc/XYZ/g' /tmp/file.txt
cp rectangle_save.poly rectangle.poly
sed -i -e 's/aa/'"$x1"'/g' rectangle.poly
sed -i -e 's/bb/'"$x2"'/g' rectangle.poly
sed -i -e 's/cc/'"$y1"'/g' rectangle.poly
sed -i -e 's/dd/'"$y2"'/g' rectangle.poly


refinement="$5"

L="$(bc <<< "scale=8;$x2 - $x1")"
# define nnx, nny, mesh size, area
nnxO="6"
nnx="$(bc <<< "scale=8;($nnxO - 1) * 2^($refinement) + 1")"
nny="$(bc <<< "scale=0;($nnx - 1)/10 + 1")"
he="$(bc <<< "scale=4;$L / ($nnx - 1)")"
area="$(bc <<< "scale=4;0.5 * $he^2")"


# read poly file and use triangle to create mesh
triangle -pq28Dena$area rectangle.poly
# create FEM file here for compatibility of Guermond's code
triangle_2_fem rectangle.1

# clean up, note this was easiest way for now without doing a better way to write
# temporary files etc

mv rectangle.1.FEM rectangle.FEM
mv rectangle.1.txt rectangle.txt
# don't do this if want to use showme
rm rectangle.1.*
rm rectangle.poly

echo "this is number of nodes in x $nnx"
echo "this is number of nodes in y $nny"
echo "this is mesh size            $he"
echo "this is area                 $area"
