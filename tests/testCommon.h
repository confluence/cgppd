#ifndef TEST_H_
#define TEST_H_

#include <set>

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/SourceLine.h>
#include <cppunit/TestAssert.h>

#include <vector3f.h>
#include <Potential.h>

void check_vector3f_equal(Vector3f expected, Vector3f actual, CppUnit::SourceLine sourceLine);
#define ASSERT_VECTOR3F_EQUALS(expected, actual) check_vector3f_equal(expected, actual, CPPUNIT_SOURCELINE())

void check_potential_equal(double * expected, Potential actual, CppUnit::SourceLine sourceLine);
#define ASSERT_POTENTIAL_EQUALS(expected, actual) check_potential_equal(expected, actual, CPPUNIT_SOURCELINE())

template <typename Type>
string join(const string sep, set<Type> values)
{
    string s;

    for (typename set<Type>::iterator i = values.begin(); i != values.end(); i++) {
        s += sep + string(*i);
    }

    return s;
}

template <>
string join(const string sep, set<int> values);


template <typename Type>
void check_set_equal(set<Type> expected, set<Type> actual, CppUnit::SourceLine sourceLine)
{
    if (expected != actual) {
        ::CppUnit::Asserter::failNotEqual(join(" ", expected), join(" ", actual), sourceLine);
    }
    
    return;
}
#define ASSERT_SET_EQUALS(expected, actual) check_set_equal(expected, actual, CPPUNIT_SOURCELINE())

#endif // TEST_H_
