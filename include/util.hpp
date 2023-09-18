// Utility string and log functions.
#include <string>

char *strcpyTrim(std::string *s, std::string *ct);
void logAndExit(FILE *logfp, std::string *outstr);
void logAndPrint(FILE *logfp, std::string *outstr);
void logOnly(FILE *logfp, std::string *outstr);
bool namesMatch(const std::string &n1, const std::string &n2);

