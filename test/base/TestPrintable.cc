#include "catch.hpp"

#include "base/printable.h"

#include <sstream>
#include <vector>
#include <list>
#include <set>
#include <map>

using namespace std;
using namespace ant;

void dotest();

TEST_CASE("Printable", "[base]") {
    dotest();
}

struct A : printable_traits {
    virtual std::ostream& Print( std::ostream& stream ) const {
        return stream << "A";
    }
        ;
};

void dotest() {
    A a;
    {
        stringstream ss;
        ss << a;
        CHECK(ss.str() == "A");
    }

    {
        vector<A> c_a{};
        stringstream ss;
        ss << c_a;
        CHECK(ss.str() == "[]");
    }

    {
        vector<A> c_a{a};
        stringstream ss;
        ss << c_a;
        CHECK(ss.str() == "[A]");
    }

    {
        vector<A> c_a{a,a};
        stringstream ss;
        ss << c_a;
        CHECK(ss.str() == "[A;A]");
    }

    {
        list<A> c_a{};
        stringstream ss;
        ss << c_a;
        CHECK(ss.str() == "[]");
    }

    {
        list<A> c_a{a};
        stringstream ss;
        ss << c_a;
        CHECK(ss.str() == "[A]");
    }

    {
        list<A> c_a{a,a};
        stringstream ss;
        ss << c_a;
        CHECK(ss.str() == "[A;A]");
    }

    {
        map<unsigned, A> c_a{};
        stringstream ss;
        ss << c_a;
        CHECK(ss.str() == "[]");
    }

    {
        map<unsigned, A> c_a;
        c_a[0] = a;
        stringstream ss;
        ss << c_a;
        CHECK(ss.str() == "[0=A]");
    }

    {
        map<unsigned, A> c_a;
        c_a[0] = a;
        c_a[1] = a;
        stringstream ss;
        ss << c_a;
        CHECK(ss.str() == "[0=A;1=A]");
    }

    {
        set<unsigned> c_a{};
        stringstream ss;
        ss << c_a;
        CHECK(ss.str() == "[]");
    }

    {
        set<unsigned> c_a{1};
        stringstream ss;
        ss << c_a;
        CHECK(ss.str() == "[1]");
    }

    {
        set<unsigned> c_a{1,2};
        stringstream ss;
        ss << c_a;
        CHECK(ss.str() == "[1;2]");
    }

    {
        set<unsigned> c_a{1,1};
        stringstream ss;
        ss << c_a;
        CHECK(ss.str() == "[1]");
    }

}
