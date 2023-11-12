cd "/mnt/c/Programming/Git Projects/tsunami_lab"

rm -R "build"
scons -Q debug=0
cd "build"

cd "/mnt/c/Programming/Git Projects/tsunami_lab"

rm -R "html"
rm -R "latex"
doxygen

cd "/mnt/c/Programming/Git Projects/tsunami_lab/build"

./tests
./tsunami_lab 100 FWAVE BATHYMETRY OUTFLOW OUTFLOW