cd "/mnt/c/Programming/Git Projects/tsunami_lab"

rm -R "build"
scons -Q debug=0
cd "build"

./tests
./tsunami_lab 10 FWAVE

cd "/mnt/c/Programming/Git Projects/tsunami_lab"

rm -R "html"
rm -R "latex"
doxygen