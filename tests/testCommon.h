#ifndef TEST_H_
#define TEST_H_

#include <set>
#include <vector>

#define ELPP_NO_DEFAULT_LOG_FILE 1
#include "easylogging++.h"
#include "catch.hpp"

using namespace std;

vector<int> make_range(int begin, int end);
set<int> to_set(vector<int> v);

#endif // TEST_H_
