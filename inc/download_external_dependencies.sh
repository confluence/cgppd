#/!/bin/bash

wget https://github.com/easylogging/easyloggingpp/releases/download/v9.80/easyloggingpp_v9.80.tar.gz &&
tar -xzf easyloggingpp_v9.80.tar.gz &&
rm easyloggingpp_v9.80.tar.gz &&
chmod -x "easylogging++.h" &&
# I know this is horrifying, but there's a macro name clash with catch.hpp
sed -ri 's/\bCHECK\b/ELPP_CHECK/g' "easylogging++.h"

wget https://raw.githubusercontent.com/philsquared/Catch/develop/single_include/catch.hpp
