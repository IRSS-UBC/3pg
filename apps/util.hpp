// Utility string and log functions.
#include <string>
using namespace std;

char *strcpyTrim(char *s, char *ct);
bool namesMatch(const std::string &n1, const std::string &n2);

class Logger
{
private:
	string logName;
public:
	Logger(const string& filename);
	~Logger();
	string GetCurrentDate();
	string GetCurrentTime();
	void Log(const string& logMsg);
};
