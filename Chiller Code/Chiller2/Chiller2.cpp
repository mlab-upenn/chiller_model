// Chiller2.cpp : Defines the entry point for the console application.
//

#include "StdAfx.h"
#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "Headers\Protos.h"
#include "Headers\Properties.h"
#include "Headers\Components.h"
#include "Headers\NumRecipes.h"
#include "Headers\ErrorLog.h"

//////////////////////////
#define	SS1		2000
#define SS2		3500
#define SS3		5500
#define	SS4		7500
#define	SS5		9500
#define	SS6		11000
#define	SS7		13000
#define	SS8		15000
#define	SS9		16500
#define	SS10	19000
#define	SS11	20500
#define	SS12	22500
#define	SS13	24500
#define	SS14	26000
#define SS15	28000
#define	SS16	29500
#define	SS17	31500
#define	SS18	33500
#define	SS19	35500
#define	SS20	37500
#define	SS21	39500
#define	SS22	41000
#define	SS23	43000
#define	SS24	44500
#define	SS25	46000
#define	SS26	48000
#define	SS27	49500
////////////////////////

ErrorLog	errorlog;

//Global System instance
VapCompCentLiqChiller*	chiller = NULL;

#define	RUNSIM	200.0

bool	TisSS(int t);
int	main(int argc, char *argv[])
{

	if(chiller==NULL)
	{
		chiller	= new VapCompCentLiqChiller();
		if(chiller==NULL)
		{
			cout << "Cannot create engine.";
			errorlog.Add("Main::VapCompCentLiqChiller ","Object constructor failed.");
			return 1;
		}
		if(!chiller->fnDefineGeometry())
		{
			errorlog.Add("Main::Chiller->fnDefineGeometry ","System construction failed.");
			return 1;
		}
	}

	if(errorlog.IsError())		errorlog.ClearError();
	if(!chiller->fnLoadState(MINIMAL))
	{
		errorlog.Add("Main::Chiller->fnLoadState ","System initialization failed.");
		return 1;
	}

	chiller->VINPUTS(1,16);	chiller->VINPUTS(2,29);	chiller->VINPUTS(3,10);	chiller->VINPUTS(4,13);	chiller->VINPUTS(5,17);

	ifstream	infile("IOFiles\\FFBCs.txt");
	ofstream	outfile("IOFiles\\FF_Output.txt");
	double	Tewin, Tcwin, Tewoset, mewat, mcwat;
	clock_t	start, stepstart, stepend;
	start = clock();
	double RTF=0.0, ORTF=0.0;
	char errstr[512];
	while(chiller->dTIME<=RUNSIM)
	{
		infile >> Tewin >> Tcwin >> Tewoset >> mewat >> mcwat;
		chiller->VINPUTS(1,Tewin);	chiller->VINPUTS(2,Tcwin);	chiller->VINPUTS(3,Tewoset);
		chiller->VINPUTS(4,mewat);	chiller->VINPUTS(5,mcwat);
		for(int i=1;i<=10;i++)
		{
			printf("Step %d\tof %d at %6.4fs per step\tOverall:%f%c",(int)chiller->dTIME,(int)RUNSIM,RTF,ORTF,13);
			stepstart = clock();
			if(!chiller->fnAdvance1sec(EXPEPC))
			{
				sprintf(errstr,"Chiller execution failed at time %f",chiller->dTIME);
				errorlog.Add("Main::Chiller->fnAdvance1sec ",errstr);
				infile.close();
				outfile.close();
				return 1;
			}
			stepend = clock();
			RTF = (stepend-stepstart)/(double)CLOCKS_PER_SEC;
			ORTF = (stepend-start)/(chiller->dTIME*(double)CLOCKS_PER_SEC);
			outfile << chiller->dTIME << "\t";
			for(int j=1;j<=chiller->iN_OUTPUTS;j++)
				outfile << chiller->VOUTPUTS(j)	<< "\t";
			outfile << RTF << endl;
		}
		if(!chiller->fnSaveState())
		{
			sprintf(errstr,"Chiller state could not be saved at time %f",chiller->dTIME);
			errorlog.Add("Main::Chiller->fnSaveState ",errstr);
		}
		if(TisSS((int)chiller->dTIME))
		{
			char fname[512];
			sprintf(fname,"SystemState_%d.txt",(int)chiller->dTIME);
			if(!chiller->fnSaveState(fname))
			{
				sprintf(errstr,"Chiller state could not be saved to file %s",fname);
				errorlog.Add("Main::Chiller->fnSaveState ",errstr);
			}
		}
	}


	infile.close();
	outfile.close();
	return 0;
}

bool	TisSS(int t)
{
	return (t==SS1||t==SS2||t==SS3||t==SS4||t==SS5||t==SS6||t==SS7||t==SS8||t==SS9||
		t==SS10||t==SS11||t==SS12||t==SS13||t==SS14||t==SS15||t==SS16||t==SS17||t==SS18||
		t==SS19||t==SS20||t==SS21||t==SS22||t==SS23||t==SS24||t==SS25||t==SS26||t==SS27);

}


