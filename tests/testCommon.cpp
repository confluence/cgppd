#include <testCommon.h>

void check_vector3f_equal(Vector3f expected, Vector3f actual, CppUnit::SourceLine sourceLine)
{
    float eps(0.00001);
    if (fabs(expected.x - actual.x) < eps && fabs(expected.y - actual.y) < eps && fabs(expected.z - actual.z) < eps)
    {
        return;
    }

    char expected_char[256];
    sprintf(expected_char, "vector3f(%f, %f, %f)", expected.x, expected.y, expected.z);
    string expected_str(expected_char);

    char actual_char[256];
    sprintf(actual_char, "vector3f(%f, %f, %f)", actual.x, actual.y, actual.z);
    string actual_str(actual_char);

    ::CppUnit::Asserter::failNotEqual(expected_str, actual_str, sourceLine);
}

void check_potential_equal(double * expected, Potential actual, CppUnit::SourceLine sourceLine)
{
    float eps(0.00001);

    if (fabs(expected[0] - actual.total_LJ()) < eps &&
        fabs(expected[1] - actual.total_DH()) < eps &&
#if FLEXIBLE_LINKS
        fabs(expected[2] - actual.total_bond()) < eps &&
        fabs(expected[3] - actual.total_angle()) < eps &&
        fabs(expected[4] - actual.total_torsion()) < eps &&
#endif // FLEXIBLE_LINKS
        fabs(expected[5] - actual.total()) < eps)
    {
        return;
    }

    char expected_char[256];
#if FLEXIBLE_LINKS
    sprintf(expected_char, "LJ: %f DH: %f bond: %f angle: %f torsion: %f TOTAL: %f", expected[0], expected[1], expected[2], expected[3], expected[4], expected[5]);
#else
    sprintf(expected_char, "LJ: %f DH: %f TOTAL: %f", expected[0], expected[1], expected[5]);
#endif // FLEXIBLE_LINKS
    string expected_str(expected_char);

    char actual_char[256];
    actual.to_string(actual_char);
    string actual_str(actual_char);

    ::CppUnit::Asserter::failNotEqual(expected_str, actual_str, sourceLine);
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
