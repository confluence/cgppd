#include <testCommon.h>

vector<int> make_range(int begin, int end)
{
    vector<int> v;
    for (int i = begin; i <= end; i++) {
        v.push_back(i);
    }
    return v;
}

set<int> to_set(vector<int> v)
{
    set<int> s(v.begin(), v.end());
    return s;
}
