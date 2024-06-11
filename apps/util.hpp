// Utility string and log functions.
#include <string>

char *strcpyTrim(char *s, char *ct);
void logAndExit(FILE *logfp, char *outstr);
void logAndPrint(FILE *logfp, char *outstr);
void logOnly(FILE *logfp, char *outstr);
bool namesMatch(const std::string &n1, const std::string &n2);

