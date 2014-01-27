#include "StdAfx.h"
#include <time.h>
#include "Headers\errorlog.h"

ErrorLog errorLog;

ErrorLog::ErrorLog()
{
	isError=0;

	int cnt = 1;
	FILE* fp;
	do {
		sprintf(fn,"error%d.log",cnt++);
		fp = fopen(fn,"w");
	} while(fp==NULL && cnt<100);
	if(fp) {
		char da[32],ti[32];
		_strdate(da);
		_strtime(ti);
		fprintf(fp,"CentChillerSim: Starting error log %s %s.\n\n",da,ti);
		fclose(fp);
	}
}

ErrorLog::~ErrorLog()
{
}

void
ErrorLog::ClearError(const char* msg)
{
	isError = 0;
	//myErrorNumber = 0;
	FILE* fp = fopen(fn,"a");
	if(fp) {
		fprintf(fp,"Error Cleared. %s\n\n",msg);
		fclose(fp);
	}
}

void
ErrorLog::AddWarning(const char* CLASS, const char* METHOD, const char* msg)
{
	//if(fp) fprintf(fp,"Warning! %s::%s >>> %s\n",CLASS,METHOD,msg);
}

void
ErrorLog::Add(const char* func, const char* msg)
{
	isError=1;
	FILE* fp = fopen(fn,"a");
	if(fp) {
		fprintf(fp,"Fatal! %s >>> %s\n",func,msg);
		fclose(fp);
	}
}

void
ErrorLog::Add(const char* func, const char* msg, double values)
{
	isError=1;
	FILE* fp = fopen(fn,"a");
	if(fp) {
		fprintf(fp,"Fatal! %s >>> %s \t %f\n",func,msg,values);
		fclose(fp);
	}
}
