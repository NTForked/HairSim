
#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include <iostream>

/** Converts from a string to the given type */
template<class T>
inline void fromString(T& t, const std::string& str)
{
    std::stringstream ss(str);
    ss >> t;
}

/** Converts the given input to a string */
template<class T>
inline std::string toString(const T& t)
{
    std::stringstream ss;
    ss << t;
    return ss.str();
}

#endif
