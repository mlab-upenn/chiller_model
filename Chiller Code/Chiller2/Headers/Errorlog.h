#ifndef _ERRORLOG
#define _ERRORLOG

#include <stdio.h>

class ErrorLog {
public:
	ErrorLog();
	~ErrorLog();
	void Add(const char* func, const char* msg = "");
	void Add(const char* func, const char* msg, double values);
	void AddWarning(const char* CLASS, const char* METHOD, const char* msg = "");
	int IsError() {return isError;}
	void ClearError(const char* msg = "");
private:
	int isError;
	char fn[256];
};

#endif