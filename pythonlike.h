#ifndef __pythonlike_h__
#define __pythonlike_h__

/*
Author: Yotam Gingold <yotam (strudel) yotamgingold.com>
License: Public Domain [CC0](http://creativecommons.org/publicdomain/zero/1.0/)
Adapted from my old "stl.h"
On GitHub as a gist: https://gist.github.com/yig/32fe51874f3911d1c612
*/

#include <map>
#include <vector>
#include <utility> // pair
#include <algorithm> // find
#include <string>
#include <sstream>
#include <numeric>
#include <sys/stat.h> // stat() used in os_path_exists()

namespace pythonlike
{

// Behaves like the Python function of the same name.
template< typename T >
inline bool in( const std::vector< T >& vec, const T& item )
{
    return vec.end() != std::find( vec.begin(), vec.end(), item );
}

// Behaves like the Python function of the same name.
template< typename T, typename U >
inline void keys( const std::map< T, U >& map, std::vector< T >& result )
{
    result.clear();
    result.reserve( map.size() );
    
    for( typename std::map< T, U >::const_iterator it = map.begin(); it != map.end(); ++it )
    {
        result.push_back( it->first );
    }
}
template< typename T, typename U >
inline std::vector< T > keys( const std::map< T, U >& map )
{
    std::vector< T > result;
    keys( map, result );
    return result;
}

// Behaves like the Python function of the same name.
template< typename T, typename U >
inline void values( const std::map< T, U >& map, std::vector< T >& result )
{
    result.clear();
    result.reserve( map.size() );
    
    for( typename std::map< T, U >::const_iterator it = map.begin(); it != map.end(); ++it )
    {
        result.push_back( it->second );
    }
}
template< typename T, typename U >
inline std::vector< T > values( const std::map< T, U >& map )
{
    std::vector< T > result;
    values( map, result );
    return result;
}

// Behaves like the Python function of the same name.
template< typename T, typename U >
inline void items( const std::map< T, U >& map, std::vector< std::pair< T, U > >& result )
{
    result.clear();
    result.reserve( map.size() );
    
    for( typename std::map< T, U >::const_iterator it = map.begin(); it != map.end(); ++it )
    {
        result.push_back( *it );
    }
}
template< typename T, typename U >
inline void items( const std::map< T, U >& map, std::vector< T >& result_keys, std::vector< T >& result_vals )
{
    result_keys.clear();
    result_vals.clear();
    result_keys.reserve( map.size() );
    result_vals.reserve( map.size() );
    
    for( typename std::map< T, U >::const_iterator it = map.begin(); it != map.end(); ++it )
    {
        result_keys.push_back( it->first );
        result_keys.push_back( it->second );
    }
}
template< typename T, typename U >
inline std::vector< std::pair< T, U > > items( const std::map< T, U >& map )
{
    std::vector< std::pair< T, U > > result;
    result.reserve( map.size() );
    
    for( typename std::map< T, U >::const_iterator it = map.begin(); it != map.end(); ++it )
    {
        result.push_back( *it );
    }
    
    return result;
}

// Behaves like the python zip() function.
template< typename T, typename U >
inline std::vector< std::pair< T, U > > zip( const std::vector< T >& first, const std::vector< U >& second )
{
    const unsigned int min_size = std::min( first.size(), second.size() );
    std::vector< std::pair< T, U > > result;
    result.reserve( min_size );
    for( unsigned int i = 0; i < min_size; ++i )
    {
        result.push_back( std::make_pair( first.at(i), second.at(i) ) );
    }
    return result;
}

// The inverse of zip, with output parameters.
// unzip( zip( a, b ), a, b ) leaves a and b unchanged.
template< typename T, typename U >
inline void unzip( const std::vector< std::pair< T, U > >& both, std::vector< T >& first, std::vector< U >& second )
{
    first.clear();
    second.clear();
    first.reserve( both.size() );
    second.reserve( both.size() );
    
    for( unsigned int i = 0; i < both.size(); ++i )
    {
        first.push_back( both.at(i).first );
        second.push_back( both.at(i).second );
    }
}
// The inverse of zip, returning a pair of vectors instead of a vector of pairs.
// std::make_pair( a,b ) == unzip( zip( a, b ) ).
template< typename T, typename U >
inline std::pair< std::vector< T >, std::vector< U > > unzip( const std::vector< std::pair< T, U > >& both )
{
    std::pair< std::vector< T >, std::vector< U > > result;
    unzip( both, result.first, result.second );
    return result;
}

// Behaves like the python os.path.split() function.
inline std::pair< std::string, std::string > os_path_split( const std::string& path )
{
    const std::string::size_type split = path.find_last_of( "/" );
    if( split == std::string::npos )
        return std::make_pair( std::string(), path );
    else
    {
        std::string::size_type split_start = split;
        // Remove multiple trailing slashes.
        while( split_start > 0 && path[ split_start-1 ] == '/' ) split_start -= 1;
        // Don't remove the leading slash.
        if( split_start == 0 ) split_start = 1;
        return std::make_pair( path.substr( 0, split_start ), path.substr( split+1 ) );
    }
}

// Behaves like the python os.path.splitext() function.
inline std::pair< std::string, std::string > os_path_splitext( const std::string& path )
{
    const std::string::size_type split_dot = path.find_last_of( "." );
    const std::string::size_type split_slash = path.find_last_of( "/" );
    if( split_dot != std::string::npos && (split_slash == std::string::npos || split_slash < split_dot) )
        return std::make_pair( path.substr( 0, split_dot ), path.substr( split_dot ) );
    else
        return std::make_pair( path, std::string() );
}

// Behaves like the python os.path.exists() function.
inline bool os_path_exists( const std::string& path )
{
    struct stat buffer;
    // stat() returns zero for success (file exists).
    return 0 == stat( path.c_str(), &buffer );
}


// Can be used with std::transform() to transform a range using an std::map
template< typename T, typename U >
struct mapper
{
    mapper( std::map< T, U >& amap ) : map( amap ) {}
    
    U& operator()( const T& key ) { return map[ key ]; }
    
private:
    std::map< T, U >& map;
};

// Converts a string to the templated type.
template< typename T >
inline T strto( const std::string& str )
// Explicitly constructing the result (result = T(); instead of result;) means that
// built-in types (int, float, etc) will default to 0, and so return that in case of error.
{ std::istringstream converter( str ); T result = T(); converter >> result; return result; }
template< typename T >
inline T strto( const std::string& str, bool& success )
// Explicitly constructing the result (result = T(); instead of result;) means that
// built-in types (int, float, etc) will default to 0, and so return that in case of error.
// Optional output parameter `success` tells the caller this explicitly.
{ std::istringstream converter( str ); T result = T(); success = bool(converter >> result); return result; }

// Converts a string to a vector of the templated type.
template< typename T >
inline std::vector<T> strtovec( const std::string& str )
{
    std::istringstream converter( str );
    std::vector<T> result;
    while( !( converter >> std::ws ).eof() )
    {
        T next;
        converter >> next;
        // Stop on first non-T word.
        // NOTE: If we want to read all T words, then change the break to continue.
        if( !converter ) break;
        result.push_back( next );
    }
    return result;
}

// Unpack an iterable to a set of variables (C++11).
// For convenience, returns an iterator to the end of the subsequence that was
// unpacked (begin + number of variable arguments).
/*
// Example:
std::vector<std::string> stuff{ "one", "two", "three" };
std::string one("a"), two("b"), three("c");
unpack( stuff.begin(), one, two, three );
std::cout << "one: " << one << std::endl;
std::cout << "two: " << two << std::endl;
std::cout << "three: " << three << std::endl;
*/ 
template< typename Iter, typename T >
Iter unpack( Iter begin, T& next ) {
	next = *begin;
	++begin;
	return begin;
}
template< typename Iter, typename T, typename... Args >
Iter unpack( Iter begin, T& next, Args&... args )
{
	next = *begin;
	++begin;
	return unpack( begin, args... );
}


// Find the sorted indices of a given vector
template <typename T>
inline std::vector<size_t> sort_indexes(const std::vector<T> &v) {

	// initialize original index locations
	std::vector<size_t> idx(v.size());
	std::iota(idx.begin(), idx.end(), 0);

	// sort indexes based on comparing values in v
	std::sort(idx.begin(), idx.end(),
	   [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

	return idx;
}

// Finds and removes an optional argument.
// Returns true if the argument was found.
// Returns false if the argument is not found or if the name is the last item in the vector
// (which could be considered an error).
// If the given name is not found, `value` stays unmodified. That way it can be
// initialized with the default value.
/*
// Example:
vector<string> args{ "foo", "bar", "--color", "red" };
string col = "blue"; // default
get_optional_parameter( args, "--color", value );
// Now `col` is "red" and `args` is { "foo", "bar" }.
*/
/*
// Protip: Collect argv[1:] into a vector named args:
vector<string> args( argv + 1, argv + argc );
*/
inline bool get_optional_parameter( std::vector< std::string >& args, const std::string& name, std::string& value )
{
    // Look for `name`.
    auto has_name = std::find( args.begin(), args.end(), name );
    
    // It's an error to find `name` as the last argument,
    // because it must be followed by the value.
    if( has_name + 1 == args.end() ) return false;
    
    // If we don't find it, return.
    if( has_name == args.end() ) return false;
    
    // Grab the following parameter as the setting_path.
    value = *( has_name + 1 );
    
    // Remove them from the arguments.
    args.erase( has_name, has_name + 2 );
    
    return true;
}
// Finds and removes an optional flag.
// Returns true if `name` is found in `args`, false otherwise.
/*
// Example:
vector<string> args{ "foo", "bar", "--clobber" };
bool clobber = get_optional_parameter( args, "--clobber" );
// Now `args` is { "foo", "bar" }.
*/
inline bool get_optional_parameter( std::vector< std::string >& args, const std::string& name )
{
    // Look for `name`.
    auto has_name = std::find( args.begin(), args.end(), name );
    
    // If we don't find it, return false.
    if( has_name == args.end() ) return false;
    
    // Remove it from the arguments.
    args.erase( has_name );
    
    return true;
}

}


/* === tester ===
// c++ -o stl_foo stl_foo.cpp

#include "/Users/yotam/Work/meshopt/git-devel/meshopt/src/common/stl.h"

#include <iostream>

void usage( const char* bin )
{
    std::cerr << "Usage: " << bin << " <split | splitext> path\n";
}

int main( int argc, char* argv[] )
{
    if( argc != 3 )
    {
        usage( argv[0] );
        return -1;
    }
    
    const std::string cmd( argv[1] );
    const std::string path( argv[2] );
    if( cmd == "split" )
    {
        std::pair< std::string, std::string > result = stl::os_path_split( path );
        std::cout << "os_path_split( '" << path << "' ): ('" << result.first << "', '" << result.second << "')\n";
    }
    else if( cmd == "splitext" )
    {
        std::pair< std::string, std::string > result = stl::os_path_splitext( path );
        std::cout << "os_path_splitext( '" << path << "' ): ('" << result.first << "','" << result.second << "')\n";
    }
    else
    {
        usage( argv[0] );
        return -1;
    }
    
    return 0;
}
*/

#endif /* __pythonlike_h__ */
