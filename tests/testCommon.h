#ifndef TEST_H_
#define TEST_H_

#include <set>
#include <vector>

#define LOGOG_USE_PREFIX 1
#include "logog/logog.hpp"

#include "catch.hpp"

using namespace std;

vector<int> make_range(int begin, int end);
set<int> to_set(vector<int> v);

#endif // TEST_H_
