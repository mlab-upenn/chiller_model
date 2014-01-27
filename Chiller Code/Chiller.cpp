// Chiller.cpp : Defines the entry point for the console application.
//

#include <math.h>
#include <mex.h>
#include "Headers\Protos.h"
#include "Headers\Properties.h"
#include "Headers\Components.h"
#include "Headers\NumRecipes.h"
#include "Headers\ErrorLog.h"

ErrorLog	errorlog;

#define	N_NODES_C	15
#define N_NODES_E	15

//	Input arguments
#define	T_IN	prhs[0]
#define	U_IN	prhs[1]

//	Output arguments
#define	Y_OUT	plhs[0]
#define X_OUT	plhs[1]

//Global System instance
VapCompCentLiqChiller*	chiller = NULL;

//USAGE: y = Chiller(t,u)
//
//	y	=	outputs
//	t	=	time
//	u	=	inputs
//
// Each call to Chiller runs the simulation for t secs, using inputs u, held constant over this time.
// The outputs of the system are passed out at the end of that time.

void	mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *t, *u, *y;
	unsigned int m,n;
	int i;

	if((nrhs!=1&&nrhs!=2)||(nrhs==1&&nlhs!=0)||(nrhs==2&&nlhs!=1))
	{
		mexErrMsgTxt("Invalid call to chiller routine.");
		return;
	}

	//Instantiate a chiller if none exists, clear existing errorlog
	if(chiller==NULL)
	{
		chiller	= new VapCompCentLiqChiller();
		if(chiller==NULL)
		{
			mexErrMsgTxt("Cannot create engine.");
			return;
		}
		if(!chiller->fnDefineGeometry())
		{
			errorlog.Add("System construction failed.");
			return;
		}
		if(errorlog.IsError())		errorlog.ClearError();
	}

	if(nrhs==1)	//Initialization/save/restore mode
	{
		m	=	mxGetM(T_IN);
		n	=	mxGetN(T_IN);
		if(!mxIsNumeric(T_IN)||mxIsComplex(T_IN)||mxIsSparse(T_IN)||!mxIsDouble(T_IN)||m!=1||n!=1)
			mexErrMsgTxt("Initialization code needs to be an integer.");

		t = mxGetPr(T_IN);
		if((int)t[0]!=0&&(int)t[0]!=1&&(int)t[0]!=2&&(int)t[0]!=3)
			mexErrMsgTxt("Invalid argument value.");

		if((int)t[0]==2)	//Save current state
		{
			if(!chiller->fnSaveState("IOFiles\\SavedState.txt"))
			{
				errorlog.Add("System state could not be saved.");
				mexErrMsgTxt("System state could not be saved.");
				return;
			}
		}
		else if((int)t[0]==3)	//Restore saved state
		{
			if(!chiller->fnLoadState("IOFiles\\SavedState.txt"))
			{
				errorlog.Add("System state could not be loaded.");
				mexErrMsgTxt("System state could not be loaded.");
				return;
			}
		}
		else					//Initialize
		{
			if(chiller->fnInitialized())
				mexWarnMsgTxt("Existing initialization is now overwritten.");
			if(!chiller->fnLoadState((int)t[0]))
			{
				errorlog.Add("System initialization failed.");
				mexErrMsgTxt("System initialization failed.");
				return;
			}
		}
	}
	else	//Execution mode
	{
		//Check if system has been initalized
		if(!chiller->fnInitialized())
			mexErrMsgTxt("Chiller not initialized.");

		//Check the first argument and make sure it is the correct format for time.
		m	=	mxGetM(T_IN);
		n	=	mxGetN(T_IN);
		if(!mxIsNumeric(T_IN)||mxIsComplex(T_IN)||mxIsSparse(T_IN)||!mxIsDouble(T_IN)||m!=1||n!=1)
			mexErrMsgTxt("t(time) is required to be a 1x1 vector.");

		//Assign the pointer to t(time);
		t	=	mxGetPr(T_IN);

		//Check the second argument and make sure it is the correct format for inputs
		m	=	mxGetM(U_IN);
		n	=	mxGetN(U_IN);
		if(!mxIsNumeric(U_IN)||mxIsComplex(U_IN)||mxIsSparse(U_IN)||!mxIsDouble(U_IN)||m!=(unsigned int)chiller->iN_INPUTS||n!=1)
		{
			char msg[128];
			sprintf(msg,"u(inputs) is required to be a %d x 1 vector.",chiller->iN_INPUTS);
			mexErrMsgTxt(msg);
		}

		//Assign the pointer to u(inputs)
		u	=	mxGetPr(U_IN);

		if(t[0]<=0)
			mexErrMsgTxt("t(time) cannot be negative.");
		else
		{
			double tt=1;
			for(int i=0;i<chiller->iN_INPUTS;i++)
				chiller->VINPUTS(i+1,u[i]);
			while(tt++ <= t[0])
				if(!chiller->fnAdvance1sec(EXPEPC))
				{
					mexErrMsgTxt("ChillerSim Error.  Consult error log file for details.");
				}
		}
		Y_OUT	=	mxCreateDoubleMatrix(chiller->iN_OUTPUTS+1,1,mxREAL);
		y		=	mxGetPr(Y_OUT);
		y[0]	=	chiller->dTIME;
		for(i=1;i<=chiller->iN_OUTPUTS;i++)
			y[i]	= chiller->VOUTPUTS(i);
// SAVE STATE COMMENTED
//		if(!chiller->fnSaveState())
//			errorlog.Add("Could not save Chiller state.");
	}
}
