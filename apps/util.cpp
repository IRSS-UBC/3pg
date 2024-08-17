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
//Logger::Logger(const string& filename)
//{
//    logName = filename;
//}
//
//void Logger::StartLog(const string& outPath)
//{
//    this->logging = true;
//    logLoc = outPath;
//    log.open(logLoc + logName, ios::trunc);
//    string currDate = GetCurrentDate();
//    string currTime = GetCurrentTime();
//    log << "-------------------\n";
//    log << "OS date: " << setw(20) << currDate << "\n";
//    log << "OS time: " << setw(20) << currTime << "\n";
//    log << "-------------------\n";
//    log.close();
//}
//
//Logger::~Logger()
//{
//    if (log.is_open())
//    {
//        log.close();
//    }
//}
//
//void Logger::Log(const std::string& logMsg)
//{
//    if (!this->logging) {
//        return;
//    }
//    log.open(logLoc + logName, std::ios::app);
//    log << logMsg + "\n";
//    log.close();
//}
//
//std::string Logger::GetCurrentDate()
//{
//    if (!this->logging) {
//        return "ERROR logger not started";
//    }
//    auto now = std::chrono::system_clock::now();
//    auto local_time = std::chrono::current_zone()->to_local(now);
//    return std::format("{:%d-%m-%Y}", local_time);
//}
//
//std::string Logger::GetCurrentTime()
//{
//    if (!this->logging) {
//        return "ERROR logger not started";
//    }
//    auto now = std::chrono::system_clock::now();
//    auto local_time = std::chrono::current_zone()->to_local(now);
//    return std::format("{:%H:%M:%S}", local_time);
//}
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
