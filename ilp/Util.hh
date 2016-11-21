/* commode -- find significant protein complex graph models
   Copyright (C) 2013, 2014  Falk Hüffner

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License along
   with this program; if not, write to the Free Software Foundation, Inc.,
   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.  */

#ifndef UTIL_HH_INCLUDED
#define UTIL_HH_INCLUDED

#include <algorithm>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace Util {

std::vector<std::string> split(const std::string& str, const std::string& delimiters = " \t\n\r");

template <typename T>
std::string join(const T& strings, const std::string& delimiter = " ") {
    std::string r;
    for (std::string s : strings) {
	if (!r.empty())
	    r += delimiter;
	r += s;
    }
    return r;
}

template<typename Container>
std::size_t intersectionSize(const Container& s1, const Container& s2) {
    std::size_t count = 0;
    auto first1 = s1.begin(), last1 = s1.end();
    auto first2 = s2.begin(), last2 = s2.end();
    while (first1 != last1 && first2 != last2) {
        if (*first1 < *first2)
            ++first1;
        else if (*first2 < *first1)
            ++first2;
        else {
            ++count;
            ++first1;
            ++first2;
        }
    }
    return count;
}

namespace SetOps {

template<typename T>
bool contains(const std::set<T>& s, const T& x) {
    return s.find(x) != s.end();
}

template<typename T>
std::set<T> singleton(const T& x) {
    std::set<T> s;
    s.insert(x);
    return s;
}

template<typename T>
T pop(std::set<T>& s) {
    T x = *s.begin();
    s.erase(s.begin());
    return x;
}

template<typename T>
std::set<T>& operator+=(std::set<T>& s, const T& x) {
    s.insert(x);
    return s;
}

template<typename T>
std::set<T>& operator-=(std::set<T>& s, const T& x) {
    s.erase(x);
    return s;
}

template<typename T>
std::set<T> operator+(const std::set<T>& s, const T& x) {
    auto r = s;
    r += x;
    return r;
}

template<typename T>
std::set<T> operator-(const std::set<T>& s, const T& x) {
    auto r = s;
    r -= x;
    return r;
}

template<typename T>
std::set<T> operator+(const std::set<T>& s1, const std::set<T>& s2) {
    std::set<T> r;
    std::set_union(s1.begin(), s1.end(), s2.begin(), s2.end(),
		   std::inserter(r, r.end()));
    return r;
}


template<typename T>
std::set<T> operator-(const std::set<T>& s1, const std::set<T>& s2) {
    std::set<T> r;
    std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(),
			std::inserter(r, r.end()));
    return r;
}

template<typename T>
std::set<T>& operator+=(std::set<T>& s, const std::set<T>& s2) {
    return s = s + s2;
}

template<typename T>
std::set<T>& operator-=(std::set<T>& s, const std::set<T>& s2) {
    return s = s - s2;
}

}

}

#endif // UTIL_HH_INCLUDED
