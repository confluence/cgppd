#ifndef TEST_H_
#define TEST_H_

#include <set>
#include <sstream>

#include "easylogging++.h"

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/SourceLine.h>
#include <cppunit/TestAssert.h>

#include <vector3f.h>
#include <Potential.h>

void check_vector3f_equal(Vector3f expected, Vector3f actual, CppUnit::SourceLine sourceLine);
#define ASSERT_VECTOR3FS_EQUAL(expected, actual) check_vector3f_equal(expected, actual, CPPUNIT_SOURCELINE())

void check_potential_equal(Potential expected, Potential actual, CppUnit::SourceLine sourceLine);
#define ASSERT_POTENTIALS_EQUAL(expected, actual) check_potential_equal(expected, actual, CPPUNIT_SOURCELINE())

template <typename Type>
string join(const string sep, Type values)
{
    string s;

    for (typename Type::iterator i = values.begin(); i != values.end(); i++) {
        s += sep + string(*i);
    }

    return s;
}

template <>
string join(const string sep, set<int> values);

template <>
string join(const string sep, vector<int> values);

template <typename Type>
void check_iterable_equal(Type expected, Type actual, CppUnit::SourceLine sourceLine)
{
    if (expected != actual) {
        ::CppUnit::Asserter::failNotEqual(join(" ", expected), join(" ", actual), sourceLine);
    }
    
    return;
}
#define ASSERT_ITERABLE_EQUALS(expected, actual) check_iterable_equal(expected, actual, CPPUNIT_SOURCELINE())

vector<int> make_range(int begin, int end);
set<int> to_set(vector<int> v);

#endif // TEST_H_
