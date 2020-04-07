# rename all files with prefix to have an extension on it

find . -name "$1*" -exec mv {} {}'.lammps.log' \;
