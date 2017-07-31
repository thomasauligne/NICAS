#!/bin/ksh
#----------------------------------------------------------------------
# Korn shell script: pack
# Author: Benjamin Menetrier
# Licensing: this code is distributed under the CeCILL-C license
# Copyright © 2017 METEO-FRANCE
#----------------------------------------------------------------------
# Make a temporary directory
cd ..
rm -fr pack
mkdir pack

# Copy LICENCE.dox, LICENSE.dox, README.dox, CHANGE_LOG.dox and .gitignore
cp -f LICENCE.dox LICENSE.dox README.dox CHANGE_LOG.dox .gitignore pack

# Copy build
mkdir pack/build
cp -f build/CMakeLists.txt pack/build

# Copy data
mkdir pack/data
links_list=`find data -name '*.ksh'`
for links in ${links_list} ; do
  cp --parents $links pack
done

# Copy doc
mkdir pack/doc
cp -f doc/Doxyfile pack/doc
cp -f doc/mainpage.h pack/doc
cp -f doc/wiki pack/doc

# Copy doc
mkdir pack/ncl
mkdir pack/ncl/script
cp -f ncl/script/*.ncl pack/ncl/script

# Copy run
mkdir pack/run
cp -f run/namelist* pack/run

# Copy script
mkdir pack/script
cp -f script/*.ksh pack/script
cp -f script/namelist.sqlite pack/script

# Copy src
mkdir pack/src
cp -f src/*.f90 pack/src
cp -f src/*.cpp pack/src
cp -f src/*.h pack/src
cp -f src/*.hpp pack/src
mkdir pack/src/external
cp -f src/external/*.f90 pack/src/external
cp -f src/external/*.c pack/src/external
cp -f src/external/*.cpp pack/src/external
cp -f src/external/*.h pack/src/external

# Copy test
mkdir pack/test
cp -f test/grid.nc pack/test
cp -f test/grid_SCRIP.nc pack/test
cp -f test/links.ksh pack/test
cp -f test/truth_param.nc pack/test

# Rename and pack everything
find pack -type f -name '*~' -delete
today=`date +%Y%m%d`
rm -fr nicas_${today}
mv pack nicas_${today}
tar -cvzf nicas_${today}.tar.gz nicas_${today}

# Clean
rm -fr nicas_${today}
