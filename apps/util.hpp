// Utility string and log functions.
#include <string>
#include <fstream>
using namespace std;

char *strcpyTrim(char *s, char *ct);
bool namesMatch(const std::string &n1, const std::string &n2);

//class Logger {
//private:
//	std::string logName;
//	std::string logLoc;
//	std::ofstream log;
//	bool logging = false;
//public:
//	Logger(const std::string& filename);
//	~Logger();
//	void StartLog(const std::string& outPath);
//	std::string GetCurrentDate();
//	std::string GetCurrentTime();
//	void Log(const std::string& logMsg);
//};
