
#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include <iostream>
#include <fstream>

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

template<typename T> 
void serializeVarHex( const T& var, std::ofstream& output_stream )
{
    assert( output_stream.good() );
    T local_var = var;
    output_stream.precision( std::numeric_limits<double>::digits10 + 2);
    output_stream.flags( std::ios_base::fixed | std::ios_base::scientific );
    output_stream.write( reinterpret_cast<char*>( &local_var ), sizeof(T) );
}


#endif
