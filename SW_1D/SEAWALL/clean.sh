echo "Cleaning up..."
mkdir tmp && mv gauge_locations.txt tmp
rm *.plt *.txt
mv tmp/* . && rm -rf tmp