// Utility string and log functions.
#include <string>
#include <fstream>
using namespace std;

char *strcpyTrim(char *s, char *ct);
bool namesMatch(const std::string &n1, const std::string &n2);

class Logger
{
private:
	string logName;
	string logLoc;
	ofstream log;
public:
	Logger(const string& filename);
	~Logger();
	void StartLog(const string& outPath);
	string GetCurrentDate();
	string GetCurrentTime();
	void Log(const string& logMsg);
};
