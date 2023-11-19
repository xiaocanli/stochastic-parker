#!/bin/bash
# fix compiling errors when using new versions of gfortran
# Copy the script to the source directory of mt_stream
# chmod +x fix_mt_stream.sh
# ./fix_mt_stream.sh

sed -i -e 's/Z'\''123'\'', Z'\''234'\'', Z'\''345'\'', Z'\''456'\''/int\(Z'\''123'\''\), int\(Z'\''234'\''\), int\(Z'\''345'\''\), int\(Z'\''456'\''\)/g' check_stream.F90
sed -i -e 's/avec = Z'\''9908b0df'\''/avec = int\(Z'\''9908b0df'\'', kind=INT32\)/g' f_jump_ahead_coeff/f_jump_coeff.F90
sed -i -e 's/ZE = Z'\''eeeeeeee'\''/ZE = int\(Z'\''eeeeeeee'\'', kind=INT32\)/g' f_jump_ahead_coeff/gf2xe.F90
sed -i -e 's/ZC = Z'\''cccccccc'\''/ZC = int\(Z'\''cccccccc'\'', kind=INT32\)/g' f_jump_ahead_coeff/gf2xe.F90
sed -i -e 's/Z8 = Z'\''88888888'\''/Z8 = int\(Z'\''88888888'\'', kind=INT32\)/g' f_jump_ahead_coeff/gf2xe.F90
sed -i -e 's/dc = Z'\''0'\''/dc = int\(Z'\''0'\'', kind=INT64\)/g' f_jump_ahead_coeff/gf2xe.F90
