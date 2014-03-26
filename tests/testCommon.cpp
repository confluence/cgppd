#include <testCommon.h>

void check_vector3f_equal(Vector3f expected, Vector3f actual, CppUnit::SourceLine sourceLine)
{
    if (expected.almost_equal(actual, 0.00001)) {
        return;
    }
    
    ostringstream s_expected;
    ostringstream s_actual;
    
    s_expected << expected;
    s_actual << actual;
    
    ::CppUnit::Asserter::failNotEqual(s_expected.str(), s_actual.str(), sourceLine);
}

void check_potential_equal(Potential expected, Potential actual, CppUnit::SourceLine sourceLine)
{
    if (expected.almost_equal(actual, 0.00001)) {
        return;
    }
    
    ostringstream s_expected;
    ostringstream s_actual;
    
    s_expected << expected;
    s_actual << actual;
    
    ::CppUnit::Asserter::failNotEqual(s_expected.str(), s_actual.str(), sourceLine);
}

template <>
string join(const string sep, set<int> values)
{
    string s;

    for (typename set<int>::iterator i = values.begin(); i != values.end(); i++) {
        char temp[11];
        sprintf(temp, "%d", *i);
        s += sep + string(temp);
    }

    return s;
}

template <>
string join(const string sep, vector<int> values)
{
    string s;

    for (typename vector<int>::iterator i = values.begin(); i != values.end(); i++) {
        char temp[11];
        sprintf(temp, "%d", *i);
        s += sep + string(temp);
    }

    return s;
}

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
