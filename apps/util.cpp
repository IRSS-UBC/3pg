/*
All source code remains the property and copyright of CSIRO. 

CSIRO accepts no responsibility for the use of 3PG(S) or of the model 3-PG in
the form supplied or as subsequently modified by third parties. CSIRO disclaims
liability for all losses, damages and costs incurred by any person as a result
of relying on this software. 
Use of this software assumes agreement to this condition of use
*/

#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <cctype>      // std::tolower
#include <algorithm>   // std::equal
#include <iostream>
#include <fstream>
#include <chrono>
#include <format>
#include <iomanip>
#include <sstream> 
#include <boost/algorithm/string/trim.hpp>
#include "util.hpp"

#ifdef WIN32
#define strncasecmp strnicmp
#endif

//--------------------------------------------------------------------------

bool ichar_equals(char a, char b)
{
    return std::tolower(static_cast<unsigned char>(a)) ==
           std::tolower(static_cast<unsigned char>(b));
}

bool iequals(const std::string& a, const std::string& b)
{
    return std::equal(a.begin(), a.end(), b.begin(), b.end(), ichar_equals);
}

bool namesMatch(const std::string& n1, const std::string& n2)
{
    // This is much easier with strings
    return n1 == n2;

  // Basically strcmp, but we compare length also so that substrings don't 
  //// match. 
  /*std::size_t l1, l2;

  l1 = n1.length();
  l2 = n2.length();
  if (l1 != l2)
    return false;
  return iequals(n1, n2);*/

}

//--------------------------------------------------------------------------
