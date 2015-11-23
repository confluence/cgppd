#/!/bin/bash

mkdir -p ~/repos inc libs

# Catch -- for unit tests

git clone git@github.com:johnwbyrd/logog.git ~/repos/logog
mkdir ~/repos/logog/build
cd ~/repos/logog/build
cmake ..
make
cd -
ln -s ~/repos/logog/include inc/logog
ln -s ~/repos/logog/build/liblogog.a libs/

# Logog -- for logging

git clone git@github.com:philsquared/Catch.git ~/repos/Catch
ln -s ~/repos/Catch/single_include/catch.hpp inc/
