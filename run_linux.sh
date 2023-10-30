cd "/mnt/c/Programming/Git Projects/tsunami_lab"

rm -R "build"
scons -Q debug=1
cd "build"

./tests
# ./tsunami_lab 10

cd "/mnt/c/Programming/Git Projects/tsunami_lab"

rm -R "html"
rm -R "latex"
doxygen