#include	"StdAfx.h"
#include	<stdlib.h>
#include	<math.h>
#include	<iostream.h>
#include	<fstream.h>
#include	"Headers\Protos.h"
#include	"Headers\NumRecipes.h"
#include	"Headers\Properties.h"
#include	"Headers\ErrorLog.h"

extern ErrorLog	errorlog;

two_phase_table R134A2PTH, R134A2PTR;
one_phase_table R134ASCTH, R134ASCTR, R134ASHTH, R134ASHTR;

two_phase_table	R122PTH, R122PTR;
one_phase_table R12SCTH, R12SCTR, R12SHTH, R12SHTR;

two_phase_table	R222PTH, R222PTR, WATERLQ, BULBTPTH;
one_phase_table R22SCTH, R22SCTR, R22SHTH, R22SHTR;

#undef	_TRAP

two_phase_table::two_phase_table()
{
	int i,j;
	rows = 10;
	cols = 10;

	table.make(rows,cols);
	splines.make(rows,cols);
	currprops.make(cols);
	PROPSRGOOD	=	false;

	if(errorlog.IsError())
	{
		errorlog.Add("two_phase_table::two_pbase_table ","Construction failed");
		return;
	}

	for(i=1;i<=rows;i++)
		for(j=1;j<=cols;j++)
		{
			table(i,j,0.0);
			splines(i,j,0.0);
		}
}

bool two_phase_table::load2ptab(char *infile)
{
	ifstream	fp(infile);
	char	msg[128];

	double temp;
	int i,j;
	Vector X,Y,Z;

	if(!fp)
	{
		errorlog.Add("two_phase_table::load2ptab ","File open failure");
		return false;
	}
	fp >> rows >> cols;

	X.make(rows);
	Y.make(rows);
	Z.make(rows);

	if(errorlog.IsError())
	{
		sprintf(msg," XYZ Memory allocation failure.  Rows = %d, Cols = %d",rows,cols);
		errorlog.Add("two_phase_table::load2ptab ",msg);
		return false;
	}

	table.make(rows,cols);
	splines.make(rows,cols);
	currprops.make(cols);
	PROPSRGOOD	=	false;

	if(errorlog.IsError())
	{
		sprintf(msg," TSC Memory allocation failure.");
		errorlog.Add("two_phase_table::load2ptab ",msg);
		return false;
	}
	//Fill the table from the datafile
	for(i=1;i<=rows;i++)
	{
		for(j=1;j<=cols;j++)
		{
			fp >> temp;
			table(i,j,temp);
		}
	}
	fp.close();


	//Build the splines table
	table.getcol(1,X);
	for(i=2;i<=cols;i++)
	{
		table.getcol(i,Y);
		spline(X,Y,rows,1e30,1e30,Z);
		splines.putcol(i,Z);
	}//for i

	return true;
}//load2ptab

Vector& two_phase_table::getprop(double x, Vector& props)
{

	//One Dimensional cubic spline interpolation function
	//x is the value in the column number xcol, which is known
	//Function returns interpolated values of all properties at that x
	//Values are stored in the currprops array.

	int i;
	Vector X,Y,Y2A;

	X.make(rows);
	Y.make(rows);
	Y2A.make(rows);
	if(errorlog.IsError())
	{
		errorlog.Add("two_phase_table::getprop ","Memory allocation failure");
		PROPSRGOOD = false;
		return props;
	}

	table.getcol(1,X);

	for(i=2;i<=cols;i++)
	{
		table.getcol(i,Y);
		splines.getcol(i,Y2A);
		currprops(i,splint(X,Y,Y2A,rows,x));
		props(i,splint(X,Y,Y2A,rows,x));
	}

	currprops(1,x);
	props(1,x);
	PROPSRGOOD = true;
	return props;
}//getprop

one_phase_table::one_phase_table()
{
	int i,j;
	rows = 10;
	cols = 10;

	table.make(rows,cols);
	splines.make(rows,cols);
	currprops.make(cols);

	if(errorlog.IsError())
	{
		PROPSRGOOD	=	false;
		errorlog.Add("one_phase_table::one_phase_table ","Memory allocation failure");
		return;
	}
	for(i=1;i<=rows;i++)
		for(j=1;j<=cols;j++)
		{
			splines(i,j,0.0);
			table(i,j,0.0);
		}
}

bool one_phase_table::load1ptab(char *infile)
{
	ifstream fp(infile);

	double temp;
	Vector X,Y,Z;

	if(!fp)
	{
		errorlog.Add("one_phase_table::load1ptab ","File open failure");
		return false;
	}
	fp >> rows >> cols;

	table.make(rows,cols);
	splines.make(rows,cols);
	currprops.make(cols);

	if(errorlog.IsError())
	{
		errorlog.Add("one_phase_table::load1ptab ","Memory allocation failure");
		return false;
	}
	//Fill the table from the datafile
	//Fix the sizes of the 1st column blocks
	int k=1;

	for(int i=1;i<=rows;i++)
	{
		for(int j=1;j<=cols;j++)
		{
			fp >> temp;
			table(i,j,temp);
		}
		if(i>1&&i<rows)
		{
			if(table(i,1)!=table(i-1,1))
				blkind[k++]=i-1;
		}
	}
	fp.close();
	blkind[k] = rows;
	nblks = k;

	//Build the splines table
	int maxblksz = blkind[1];
	for(i=2;i<=nblks;i++)
		maxblksz = (maxblksz>=blkind[i]-blkind[i-1])?maxblksz:(blkind[i]-blkind[i-1]);

	Vector X2, Y2,S2;
	int top, size;

	X2.make(maxblksz);
	Y2.make(maxblksz);
	S2.make(maxblksz);
	for(int l=3;l<=cols;l++)
	{
		for(i=1;i<=nblks;i++)
		{
			top = (i==1)?1:(blkind[i-1]+1);
			size = blkind[i] - top + 1;
			X2.reset();
			Y2.reset();
			S2.reset();
			for(int j=top,k=1;j<=blkind[i];j++,k++)
			{
				X2(k,table(j,2));
				Y2(k,table(j,l));
			}
			spline(X2,Y2,size,1e30,1e30,S2);
			for(j=top,k=1;j<=blkind[i];j++,k++)
				splines(j,l,S2(k));
		}
	}
	return true;
}//load1ptab

int one_phase_table::maxblk()
{
	int max = blkind[1];

	for(int i=2;i<=nblks;i++)
		max = (max>=blkind[i]-blkind[i-1])?max:(blkind[i]-blkind[i-1]);

	return max;
}

Vector& one_phase_table::getprop(double x1, double x2, Vector& props)
{

	//Bicubic spline properties in the first and second columns of the data table
	//In case of thermodynamic properties, these are P & h
	//In case of transport properties, these are P & T

	double y1,y2;
	Vector	X1, X2, Y1, Y2, S1, S2;
	int maxblksz = blkind[1];
	for(int i=2;i<=nblks;i++)
		maxblksz = (maxblksz>=blkind[i]-blkind[i-1])?maxblksz:(blkind[i]-blkind[i-1]);
	
	X1.make(nblks);
	X2.make(maxblksz);
	Y1.make(nblks);
	Y2.make(maxblksz);
	S1.make(nblks);
	S2.make(maxblksz);

	if(x1<table(1,1))
	{
		errorlog.Add("one_phase_table::getprop ","Pressure out of (lower) bound");
		PROPSRGOOD	=	false;
		return props;
	}
	if(x1>table(rows,1))
	{
		errorlog.Add("one_phase_table::getprop ","Pressure out of (upper) bound");
		PROPSRGOOD	=	false;
		return props;
	}
	if(x2<table(1,2))
	{
		errorlog.Add("one_phase_table::getprop ","Temperature/enthalpy out of (lower) bound");
		PROPSRGOOD	=	false;
		return props;
	}
	if(x2>table(rows,2))
	{
		errorlog.Add("one_phase_table""getprop ","Temperature/enthalpy out of (upper) bound");
		PROPSRGOOD	=	false;
		return props;
	}

	int top, size,size1;
			
	for(int l=3;l<=cols;l++)
	{
		Y1.reset();
		S1.reset();
		X1.reset();
		size1=0;
		for(i=1;i<=nblks;i++)
		{
			top = (i==1)?1:(blkind[i-1]+1);
			size = blkind[i] - top + 1;
			X2.reset();
			Y2.reset();
			S2.reset();
			if(x2<table(top,2)||x2>table(blkind[i],2))
				continue;
			for(int j=top,k=1;j<=blkind[i];j++,k++)
			{
				X2(k,table(j,2));
				Y2(k,table(j,l));
				S2(k,splines(j,l));
			}
			size1++;
			y2 = splint(X2,Y2,S2,size,x2);
			Y1(size1,y2);
			X1(size1,table(blkind[i],1));
		}//for i
		spline(X1,Y1,size1,1e30,1e30,S1);
		y1 = splint(X1,Y1,S1,size1,x1);
		currprops(l,y1);
		props(l,y1);
	}

	currprops(1,x1);
	currprops(2,x2);
	props(1,x1);
	props(2,x2);
	PROPSRGOOD	=	true;
	return props;
}//getprop

watstate::watstate()	{	UPDATE= false;	}
watstate::~watstate()	{}

bool watstate::setstate(double Temp)
{
	Vector props;

	if(!props.make(5))	return false;

	T = Temp;
	WATERLQ.getprop(T,props);
	if(!WATERLQ.goodprops())
	{
		errorlog.Add("watstate::setstate ","water properties could not be set");
		return false;
	}
	if((k	= props(2))==INF)	errorlog.Add("watstate::setstate ","conductivity not computed correctly");
	if((rho = props(3))==INF)	errorlog.Add("watstate::setstate ","density not computed correctly");
	if((mu	= props(4))==INF)	errorlog.Add("watstate::setstate ","viscosity not computed correctly");
	if((Cp	= props(5))==INF)	errorlog.Add("watstate::setstate ","specific heat not computed correctly");

	if(errorlog.IsError())
	{
		errorlog.Add("watstate::setstate ","Some/all transport properties of water not computed correctly");
		return false;
	}

	UPDATE	= true;
	return true;
}
double watstate::getT()			{	return T;	}
double watstate::getk()			{	return k;	}
double watstate::getrho()		{	return rho;	}
double watstate::getmu()		{	return mu;	}
double watstate::getCp()		{	return Cp;	}
bool   watstate::getupdate()	{	return UPDATE;	}

/////////////////    R  134A //////////////////////////////////////////////////////////////////////

int		r134astate::setphase(double enth, double hl, double hv)
{
	double	temp;
	int		temp2;

	temp2 = int(enth*1000.0);
	temp = (int(enth*1000.0))/1000.0;

	if(enth>hl&&enth<hv)	//two phase
	{
		x = (enth-hl)/(hv-hl);
		return (phase = TWOPHASE);
	}
	else if(enth<=hl)	//subcooled
	{
		x=-1.0;
		return (phase = SUBCOOLED);
	}
	else						//superheated
	{
		x=1.1;
		return (phase = SUPERHEATED);
	}
}//setphase and quality

r134astate::r134astate()	{	UPDATE = false;	}
r134astate::~r134astate() {}


bool r134astate::setstate(double Pres, double enth)
{
	//sets all the properties in the state for given pressure and enthalpy
	Vector props;
	P=Pres;
	h=enth;

	if(!props.make(12))	return false;

	//First determine all saturation properties at given pressure
	props = R134A2PTH.getprop(P,props);
	if(!R134A2PTH.goodprops())	return false;
	if((Ts = props(2))==INF) return false;
	if((hl = props(3))==INF) return false;
	if((hv = props(4))==INF) return false;
	if((vl = props(5))==INF) return false;
	if((vv = props(6))==INF) return false;
	rhol = 1/vl;
	rhov = 1/vv;
	if((sl = props(7))==INF) return false;
	if((sv = props(8))==INF) return false;
	if((ul = props(9))==INF) return false;
	if((uv = props(10))==INF) return false;
	props = R134A2PTR.getprop(P,props);
	if(!R134A2PTR.goodprops())	return false;
	if((kl = props(3))==INF) return false;
	if((kv = props(4))==INF) return false;
	if((mul= props(5))==INF) return false;
	if((muv= props(6))==INF) return false;
	if((Cpl= props(7))==INF) return false;
	if((Cpv= props(8))==INF) return false;

	{
		double dP = 1e-3;
		double hvPdP, hvP_dP, hlPdP, hlP_dP;
		double rvPdP, rvP_dP, rlPdP, rlP_dP;
		props = R134A2PTH.getprop(P+dP,props);
		if(!R134A2PTH.goodprops())	return false;
		hvPdP = props(4);
		hlPdP = props(3);
		rvPdP = 1.0/props(6);
		rlPdP = 1.0/props(5);
		props = R134A2PTH.getprop(P-dP,props);
		if(!R134A2PTH.goodprops())	return false;
		hvP_dP = props(4);
		hlP_dP = props(3);
		rvP_dP = 1.0/props(6);
		rlP_dP = 1.0/props(5);

		dhvdP = (hvPdP - hvP_dP)/(2*dP);
		dhldP = (hlPdP - hlP_dP)/(2*dP);
		drvdP = (rvPdP - rvP_dP)/(2*dP);
		drldP = (rlPdP - rlP_dP)/(2*dP);
	}

	phase = setphase(h,hl,hv);

	switch(phase)
	{
	case SUBCOOLED:			//subcooled
		if((hl-h)>50)
		{
			errorlog.Add("r134astate::setstate ","Enthalpy is out of (lower) bound");
			return false;
		}
		//...compressible...//
		props = R134ASCTH.getprop(P,h,props);
		if(!R134ASCTH.goodprops())	return false;
		if((T	= props(3))==INF)	return false;
		if((v	= props(4))==INF)	return false;
		rho	=	1/v;
		if((s	= props(5))==INF)	return false;
		if((u	= props(6))==INF)	return false;
		if((drdP= props(7))==INF)	return false;
		if((drdh= props(8))==INF)	return false;

		props = R134ASCTR.getprop(P,T,props);
		if(!R134ASCTR.goodprops())		return false;
		if((k	= props(3))==INF) return false;
		if((mu	= props(4))==INF) return false;
		if((Cp	= props(5))==INF) return false;

		break;
	case TWOPHASE:			//two-phase
		T	= Ts;
		v	= x*vv+(1-x)*vl;
		rho	= 1/v;
		s	= x*sv+(1-x)*sl;
		u	= x*uv+(1-x)*ul;
		{
			double x1, x2, dh=1e-3,dP=1e-3;
			double rho1, rho2;

			x1 = (h+dh-hl)/(hv-hl);
			x2 = (h-dh-hl)/(hv-hl);
			rho1 = 1/(x1*vv+(1-x1)*vl);
			rho2 = 1/(x2*vv+(1-x2)*vl);
			drdh = (rho1-rho2)/(2*dh);

			props = R134A2PTH.getprop(P+dP,props);
			if(!R134A2PTH.goodprops())		return false;
			x1 = (h-props(3))/(props(4)-props(3));
			rho1 = 1/(x1*props(6)+(1-x1)*props(5));
			props = R134A2PTH.getprop(P-dP,props);
			if(!R134A2PTH.goodprops())		return false;
			x2 = (h-props(3))/(props(4)-props(3));
			rho2 = 1/(x2*props(6)+(1-x2)*props(5));
			drdP = (rho1-rho2)/(2*dP);
		}
		k	= x*kv+(1-x)*kl;
		mu	= x*muv+(1-x)*mul;
		Cp	= x*Cpv+(1-x)*Cpl;
		break;
	case SUPERHEATED:			//superheated
		props = R134ASHTH.getprop(P,h,props);
		if(!R134ASHTH.goodprops())		return false;
		if((T	= props(3))==INF) return false;
		if((v	= props(4))==INF) return false;
		rho	= 1/v;
		if((s	= props(5))==INF) return false;
		if((u	= props(6))==INF) return false;
		if((drdP= props(7))==INF) return false;
		if((drdh= props(8))==INF) return false;
		props = R134ASHTR.getprop(P,T,props);
		if(!R134ASHTR.goodprops())		return false;
		if((k	= props(3))==INF) return false;
		if((mu	= props(4))==INF) return false;
		if((Cp	= props(5))==INF) return false;
		break;
	default:
		break;
	}//switch
	UPDATE = true;
	return true;
}

bool r134astate::setstateP(double Press)
{
	Vector props(12);
	P=Press;

	//Determine all saturation properties at given pressure
	props = R134A2PTH.getprop(P,props);
	if(!R134A2PTH.goodprops())	return false;
	if((Ts = props(2))==INF) return false;
	if((hl = props(3))==INF) return false;
	if((hv = props(4))==INF) return false;
	if((vl = props(5))==INF) return false;
	if((vv = props(6))==INF) return false;
	rhol = 1/vl;
	rhov = 1/vv;
	if((sl = props(7))==INF) return false;
	if((sv = props(8))==INF) return false;
	if((ul = props(9))==INF) return false;
	if((uv = props(10))==INF) return false;
	props = R134A2PTR.getprop(P,props);
	if(!R134A2PTR.goodprops())	return false;
	if((kl = props(3))==INF) return false;
	if((kv = props(4))==INF) return false;
	if((mul= props(5))==INF) return false;
	if((muv= props(6))==INF) return false;
	if((Cpl= props(7))==INF) return false;
	if((Cpv= props(8))==INF) return false;

	{
		double dP = 1e-3;
		double hvPdP, hvP_dP, hlPdP, hlP_dP;
		double rvPdP, rvP_dP, rlPdP, rlP_dP;
		props = R134A2PTH.getprop(P+dP,props);
		if(!R134A2PTH.goodprops())	return false;
		hvPdP = props(4);
		hlPdP = props(3);
		rvPdP = 1.0/props(6);
		rlPdP = 1.0/props(5);
		props = R134A2PTH.getprop(P-dP,props);
		if(!R134A2PTH.goodprops())	return false;
		hvP_dP = props(4);
		hlP_dP = props(3);
		rvP_dP = 1.0/props(6);
		rlP_dP = 1.0/props(5);

		dhvdP = (hvPdP - hvP_dP)/(2*dP);
		dhldP = (hlPdP - hlP_dP)/(2*dP);
		drvdP = (rvPdP - rvP_dP)/(2*dP);
		drldP = (rlPdP - rlP_dP)/(2*dP);
	}

	return true;
}//r134astate::setstateP

bool r134astate::setstateT(double T)
{
	double Pguess = 500.0, Ttemp,dP = 100.0, Terr=1, Ttol=1e-2;
	Vector props(12);
	bool	errsign = true, olderrsign;

	do
	{
		props = R134A2PTH.getprop(Pguess,props);
		if(!R134A2PTH.goodprops())	return false;
		Ttemp = props(2);
		Terr = T - Ttemp;
		olderrsign = errsign;
		if(Terr>0)
		{
			errsign = true;
			Pguess += dP;
		}
		else
		{
			errsign = false;
			Pguess -= dP;
		}
		if(errsign!=olderrsign)
			dP /= 10.0;
	}while(fabs(Terr)>Ttol);

	P  = Pguess;
	props = R134A2PTH.getprop(P,props);
	if(!R134A2PTH.goodprops())	return false;
	if((Ts = props(2))==INF) return false;
	if((hl = props(3))==INF) return false;
	if((hv = props(4))==INF) return false;
	if((vl = props(5))==INF) return false;
	if((vv = props(6))==INF) return false;
	rhol = 1/vl;
	rhov = 1/vv;
	if((sl = props(7))==INF) return false;
	if((sv = props(8))==INF) return false;
	if((ul = props(9))==INF) return false;
	if((uv = props(10))==INF) return false;
	props = R134A2PTR.getprop(P,props);
	if(!R134A2PTR.goodprops())	return false;
	if((kl = props(3))==INF) return false;
	if((kv = props(4))==INF) return false;
	if((mul= props(5))==INF) return false;
	if((muv= props(6))==INF) return false;
	if((Cpl= props(7))==INF) return false;
	if((Cpv= props(8))==INF) return false;

	{
		double dP = 1e-3;
		double hvPdP, hvP_dP, hlPdP, hlP_dP;
		double rvPdP, rvP_dP, rlPdP, rlP_dP;
		props = R134A2PTH.getprop(P+dP,props);
		if(!R134A2PTH.goodprops())	return false;
		hvPdP = props(4);
		hlPdP = props(3);
		rvPdP = 1.0/props(6);
		rlPdP = 1.0/props(5);
		props = R134A2PTH.getprop(P-dP,props);
		if(!R134A2PTH.goodprops())	return false;
		hvP_dP = props(4);
		hlP_dP = props(3);
		rvP_dP = 1.0/props(6);
		rlP_dP = 1.0/props(5);

		dhvdP = (hvPdP - hvP_dP)/(2*dP);
		dhldP = (hlPdP - hlP_dP)/(2*dP);
		drvdP = (rvPdP - rvP_dP)/(2*dP);
		drldP = (rlPdP - rlP_dP)/(2*dP);
	}

	T = Ts;
	return true;
}//r134astate::setstateT

bool	r134astate::setstateHS(double h, double s)
{
	r134astate	state;
	double	Pressure, Pg1, Pg2, sres1, sres2,temp, diff;
	Pg1 = 500.0;	Pg2 = 525.0;
	int iters=0;
	bool done = false;

	do
	{
		if(!state.setstate(Pg1,h))
		{
			errorlog.Add("r134astate::setstateHS ","Could not set state with given h-s values.");
			return false;
		}
		sres1	=	s - state.gets();
		if(!state.setstate(Pg2,h))
		{
			errorlog.Add("r134astate::setstateHS ","Could not set state with given h-s values.");
			return false;
		}
		sres2	=	s -	state.gets();

		diff = fabs((sres1-sres2)/s);
		if(diff>STOL)
		{
			temp	=	Pg1 - sres1*(Pg1 - Pg2)/(sres1 - sres2);
			Pg1		=	Pg2;
			Pg2		=	temp;
			iters++;
			if(iters>100)
			{
				errorlog.Add("r134astate::setstateHS ","Convergence not achieved in 100 iterations.");
				return false;
			}
		}
		else
			done = true;
	}while(!done);
	Pressure = (Pg1+Pg2)/2;
	return this->setstate(Pressure,h);
}

bool	r134astate::setstatePV(double P, double v)
{
	r134astate state;

	double hmin, hmax;
	if(!state.setstateP(P))
	{
		errorlog.Add("r134astate::setstatePV ","Could not determine two-phase properties with given pressure.");
		return false;
	}

	if(v<state.getvl())			{	hmax = state.gethl();		hmin = hmax - 20.0;		}
	else if(v>state.getvv())	{	hmin = state.gethv();		hmax = hmin + 50.0;		}
	else						{	hmin = state.gethl();		hmax = state.gethv();	}

	double residue, hmid, vmin, vmax;
	bool done = false;
	int iters = 0;
	while(!done)
	{
		//Search by bisection
		iters++;
		state.setstate(P,hmin);
		vmin = state.getv();
		state.setstate(P,hmax);
		vmax = state.getv();
		hmid = (hmin + hmax)/2.0;
		state.setstate(P,hmid);
		residue = v - state.getv();
		if(fabs(hmin-hmax)>1e-6)
		{
			if(residue>0)	
				hmin = hmid;
			else
				hmax = hmid;
		}
		else
			done = true;
	}//while !done

	return this->setstate(P,hmid);
}

bool	r134astate::setstatePT(double P, double T)
{
	r134astate state;
	double Tsat;
	if(!state.setstateP(P))
	{
		errorlog.Add("r134astate::setstatePT ","Could not determine two-phase properties with given pressure.");
		return false;
	}
	//Check for two-phase condition
	Tsat = state.getTs();
	if((T-Tsat)<1e-6)
	{
		errorlog.Add("r134astate::setstatePT ","Pressure and temperature are not independant.");
		return false;
	}
	//Not two-phase

	double hmin, hmax;
	if(T<Tsat)			{	hmax = state.gethl();	hmin = hmax - 20.0;	}
	else				{	hmin = state.gethv();	hmax = hmin + 50.0;	}

	double residue, hmid, Tmin, Tmax;
	bool done = false;
	int iters = 0;
	while(!done)
	{
		iters++;
		state.setstate(P,hmin);
		Tmin = state.getT();
		state.setstate(P,hmax);
		Tmax = state.getT();
		hmid = (hmin + hmax)/2.0;
		state.setstate(P,hmid);
		residue = T - state.getT();
		if(fabs(hmin-hmax)>1e-6)
		{
			if(residue>0)
				hmin = hmid;
			else
				hmax = hmid;
		}
		else
			done = true;
	}//while ! done
	return this->setstate(P,hmid);

}

/////////////////    R  12  //////////////////////////////////////////////////////////////////////

int		r12state::setphase(double enth, double hl, double hv)
{
	double	temp;
	int		temp2;

	temp2 = int(enth*1000.0);
	temp = (int(enth*1000.0))/1000.0;

	if(enth>hl&&enth<hv)	//two phase
	{
		x = (enth-hl)/(hv-hl);
		return (phase = TWOPHASE);
	}
	else if(enth<=hl)	//subcooled
	{
		x=-1.0;
		return (phase = SUBCOOLED);
	}
	else						//superheated
	{
		x=1.1;
		return (phase = SUPERHEATED);
	}
}//setphase and quality

r12state::r12state()	{	UPDATE = false;	}
r12state::~r12state() {}


bool r12state::setstate(double Pres, double enth)
{
	//sets all the properties in the state for given pressure and enthalpy
	Vector props;
	P=Pres;
	h=enth;

	if(!props.make(12))	return false;

	//First determine all saturation properties at given pressure
	props = R122PTH.getprop(P,props);
	if(!R122PTH.goodprops())	return false;
	if((Ts = props(2))==INF) return false;
	if((hl = props(3))==INF) return false;
	if((hv = props(4))==INF) return false;
	if((vl = props(5))==INF) return false;
	if((vv = props(6))==INF) return false;
	rhol = 1/vl;
	rhov = 1/vv;
	if((sl = props(7))==INF) return false;
	if((sv = props(8))==INF) return false;
	if((ul = props(9))==INF) return false;
	if((uv = props(10))==INF) return false;
	props = R122PTR.getprop(P,props);
	if(!R122PTR.goodprops())	return false;
	if((kl = props(3))==INF) return false;
	if((kv = props(4))==INF) return false;
	if((mul= props(5))==INF) return false;
	if((muv= props(6))==INF) return false;
	if((Cpl= props(7))==INF) return false;
	if((Cpv= props(8))==INF) return false;

	phase = setphase(h,hl,hv);

	switch(phase)
	{
	case SUBCOOLED:			//subcooled
		if((hl-h)>50)
		{
			errorlog.Add("r12state::setstate ","Enthalpy is out of (lower) bound");
			return false;
		}
		//...compressible...//
		props = R12SCTH.getprop(P,h,props);
		if(!R12SCTH.goodprops())	return false;
		if((T	= props(3))==INF)	return false;
		if((v	= props(4))==INF)	return false;
		rho	=	1/v;
		if((s	= props(5))==INF)	return false;
		if((u	= props(6))==INF)	return false;
		if((drdP= props(7))==INF)	return false;
		if((drdh= props(8))==INF)	return false;

		props = R12SCTR.getprop(P,T,props);
		if(!R12SCTR.goodprops())		return false;
		if((k	= props(3))==INF) return false;
		if((mu	= props(4))==INF) return false;
		if((Cp	= props(5))==INF) return false;

		break;
	case TWOPHASE:			//two-phase
		T	= Ts;
		v	= x*vv+(1-x)*vl;
		rho	= 1/v;
		s	= x*sv+(1-x)*sl;
		u	= x*uv+(1-x)*ul;
		{
			double x1, x2, dh=1e-3,dP=1e-3;
			double rho1, rho2;

			x1 = (h+dh-hl)/(hv-hl);
			x2 = (h-dh-hl)/(hv-hl);
			rho1 = 1/(x1*vv+(1-x1)*vl);
			rho2 = 1/(x2*vv+(1-x2)*vl);
			drdh = (rho1-rho2)/(2*dh);

			props = R122PTH.getprop(P+dP,props);
			if(!R122PTH.goodprops())		return false;
			x1 = (h-props(3))/(props(4)-props(3));
			rho1 = 1/(x1*props(6)+(1-x1)*props(5));
			props = R122PTH.getprop(P-dP,props);
			if(!R122PTH.goodprops())		return false;
			x2 = (h-props(3))/(props(4)-props(3));
			rho2 = 1/(x2*props(6)+(1-x2)*props(5));
			drdP = (rho1-rho2)/(2*dP);
		}
		k	= x*kv+(1-x)*kl;
		mu	= x*muv+(1-x)*mul;
		Cp	= x*Cpv+(1-x)*Cpl;
		break;
	case SUPERHEATED:			//superheated
		props = R12SHTH.getprop(P,h,props);
		if(!R12SHTH.goodprops())		return false;
		if((T	= props(3))==INF) return false;
		if((v	= props(4))==INF) return false;
		rho	= 1/v;
		if((s	= props(5))==INF) return false;
		if((u	= props(6))==INF) return false;
		if((drdP= props(7))==INF) return false;
		if((drdh= props(8))==INF) return false;
		props = R12SHTR.getprop(P,T,props);
		if(!R12SHTR.goodprops())		return false;
		if((k	= props(3))==INF) return false;
		if((mu	= props(4))==INF) return false;
		if((Cp	= props(5))==INF) return false;
		break;
	default:
		break;
	}//switch
	UPDATE = true;
	return true;
}

bool r12state::setstateP(double Press)
{
	Vector props(12);
	P=Press;

	//Determine all saturation properties at given pressure
	props = R122PTH.getprop(P,props);
	if(!R122PTH.goodprops())	return false;
	if((Ts = props(2))==INF) return false;
	if((hl = props(3))==INF) return false;
	if((hv = props(4))==INF) return false;
	if((vl = props(5))==INF) return false;
	if((vv = props(6))==INF) return false;
	rhol = 1/vl;
	rhov = 1/vv;
	if((sl = props(7))==INF) return false;
	if((sv = props(8))==INF) return false;
	if((ul = props(9))==INF) return false;
	if((uv = props(10))==INF) return false;
	props = R122PTR.getprop(P,props);
	if(!R122PTR.goodprops())	return false;
	if((kl = props(3))==INF) return false;
	if((kv = props(4))==INF) return false;
	if((mul= props(5))==INF) return false;
	if((muv= props(6))==INF) return false;
	if((Cpl= props(7))==INF) return false;
	if((Cpv= props(8))==INF) return false;

	return true;
}//r12state::setstateP

bool r12state::setstateT(double T)
{
	double Pguess = 500.0, Ttemp,dP = 100.0, Terr=1, Ttol=1e-2;
	Vector props(12);
	bool	errsign = true, olderrsign;

	do
	{
		props = R122PTH.getprop(Pguess,props);
		if(!R122PTH.goodprops())	return false;
		Ttemp = props(2);
		Terr = T - Ttemp;
		olderrsign = errsign;
		if(Terr>0)
		{
			errsign = true;
			Pguess += dP;
		}
		else
		{
			errsign = false;
			Pguess -= dP;
		}
		if(errsign!=olderrsign)
			dP /= 10.0;
	}while(fabs(Terr)>Ttol);

	P  = Pguess;
	props = R122PTH.getprop(P,props);
	if(!R122PTH.goodprops())	return false;
	if((Ts = props(2))==INF) return false;
	if((hl = props(3))==INF) return false;
	if((hv = props(4))==INF) return false;
	if((vl = props(5))==INF) return false;
	if((vv = props(6))==INF) return false;
	rhol = 1/vl;
	rhov = 1/vv;
	if((sl = props(7))==INF) return false;
	if((sv = props(8))==INF) return false;
	if((ul = props(9))==INF) return false;
	if((uv = props(10))==INF) return false;
	props = R122PTR.getprop(P,props);
	if(!R122PTR.goodprops())	return false;
	if((kl = props(3))==INF) return false;
	if((kv = props(4))==INF) return false;
	if((mul= props(5))==INF) return false;
	if((muv= props(6))==INF) return false;
	if((Cpl= props(7))==INF) return false;
	if((Cpv= props(8))==INF) return false;

	T = Ts;
	return true;
}//r12state::setstateT

bool	r12state::setstateHS(double h, double s)
{
	r12state	state;
	double	Pressure, Pg1, Pg2, sres1, sres2,temp, diff;
	Pg1 = 500.0;	Pg2 = 525.0;
	int iters=0;
	bool done = false;

	do
	{
		if(!state.setstate(Pg1,h))
		{
			errorlog.Add("r12state::setstateHS ","Could not set state with given h-s values.");
			return false;
		}
		sres1	=	s - state.gets();
		if(!state.setstate(Pg2,h))
		{
			errorlog.Add("r12state::setstateHS ","Could not set state with given h-s values.");
			return false;
		}
		sres2	=	s -	state.gets();

		diff = fabs((sres1-sres2)/s);
		if(diff>STOL)
		{
			temp	=	Pg1 - sres1*(Pg1 - Pg2)/(sres1 - sres2);
			Pg1		=	Pg2;
			Pg2		=	temp;
			iters++;
			if(iters>100)
			{
				errorlog.Add("r12state::setstateHS ","Convergence not achieved in 100 iterations.");
				return false;
			}
		}
		else
			done = true;
	}while(!done);
	Pressure = (Pg1+Pg2)/2;
	return this->setstate(Pressure,h);
}

bool	r12state::setstatePV(double P, double v)
{
	r12state state;

	double hmin, hmax;
	if(!state.setstateP(P))
	{
		errorlog.Add("r12state::setstatePV ","Could not determine two-phase properties with given pressure.");
		return false;
	}

	if(v<state.getvl())			{	hmax = state.gethl();		hmin = hmax - 20.0;		}
	else if(v>state.getvv())	{	hmin = state.gethv();		hmax = hmin + 50.0;		}
	else						{	hmin = state.gethl();		hmax = state.gethv();	}

	double residue, hmid, vmin, vmax;
	bool done = false;
	int iters = 0;
	while(!done)
	{
		//Search by bisection
		iters++;
		state.setstate(P,hmin);
		vmin = state.getv();
		state.setstate(P,hmax);
		vmax = state.getv();
		hmid = (hmin + hmax)/2.0;
		state.setstate(P,hmid);
		residue = v - state.getv();
		if(fabs(hmin-hmax)>1e-6)
		{
			if(residue>0)	
				hmin = hmid;
			else
				hmax = hmid;
		}
		else
			done = true;
	}//while !done

	return this->setstate(P,hmid);
}

bool	r12state::setstatePT(double P, double T)
{
	r12state state;
	double Tsat;
	if(!state.setstateP(P))
	{
		errorlog.Add("r12state::setstatePT ","Could not determine two-phase properties with given pressure.");
		return false;
	}
	//Check for two-phase condition
	Tsat = state.getTs();
	if((T-Tsat)<1e-6)
	{
		errorlog.Add("r12state::setstatePT ","Pressure and temperature are not independant.");
		return false;
	}
	//Not two-phase

	double hmin, hmax;
	if(T<Tsat)			{	hmax = state.gethl();	hmin = hmax - 20.0;	}
	else				{	hmin = state.gethv();	hmax = hmin + 50.0;	}

	double residue, hmid, Tmin, Tmax;
	bool done = false;
	int iters = 0;
	while(!done)
	{
		iters++;
		state.setstate(P,hmin);
		Tmin = state.getT();
		state.setstate(P,hmax);
		Tmax = state.getT();
		hmid = (hmin + hmax)/2.0;
		state.setstate(P,hmid);
		residue = T - state.getT();
		if(fabs(hmin-hmax)>1e-6)
		{
			if(residue>0)
				hmin = hmid;
			else
				hmax = hmid;
		}
		else
			done = true;
	}//while ! done
	return this->setstate(P,hmid);

}

/////////////////    R  22 //////////////////////////////////////////////////////////////////////

int		r22state::setphase(double enth, double hl, double hv)
{
	double	temp;
	int		temp2;

	temp2 = int(enth*1000.0);
	temp = (int(enth*1000.0))/1000.0;

	if(enth>hl&&enth<hv)	//two phase
	{
		x = (enth-hl)/(hv-hl);
		return (phase = TWOPHASE);
	}
	else if(enth<=hl)	//subcooled
	{
		x=-1.0;
		return (phase = SUBCOOLED);
	}
	else						//superheated
	{
		x=1.1;
		return (phase = SUPERHEATED);
	}
}//setphase and quality

r22state::r22state()	{	UPDATE = false;	}
r22state::~r22state() {}


bool r22state::setstate(double Pres, double enth)
{
	//sets all the properties in the state for given pressure and enthalpy
	Vector props;
	P=Pres;
	h=enth;

	if(!props.make(12))	return false;

	//First determine all saturation properties at given pressure
	props = R222PTH.getprop(P,props);
	if(!R222PTH.goodprops())	return false;
	if((Ts = props(2))==INF) return false;
	if((hl = props(3))==INF) return false;
	if((hv = props(4))==INF) return false;
	if((vl = props(5))==INF) return false;
	if((vv = props(6))==INF) return false;
	rhol = 1/vl;
	rhov = 1/vv;
	if((sl = props(7))==INF) return false;
	if((sv = props(8))==INF) return false;
	if((ul = props(9))==INF) return false;
	if((uv = props(10))==INF) return false;
	props = R222PTR.getprop(P,props);
	if(!R222PTR.goodprops())	return false;
	if((kl = props(3))==INF) return false;
	if((kv = props(4))==INF) return false;
	if((mul= props(5))==INF) return false;
	if((muv= props(6))==INF) return false;
	if((Cpl= props(7))==INF) return false;
	if((Cpv= props(8))==INF) return false;

	phase = setphase(h,hl,hv);

	switch(phase)
	{
	case SUBCOOLED:			//subcooled
		if((hl-h)>50)
		{
			errorlog.Add("r22state::setstate ","Enthalpy is out of (lower) bound");
			return false;
		}
		//...compressible...//
		props = R22SCTH.getprop(P,h,props);
		if(!R22SCTH.goodprops())	return false;
		if((T	= props(3))==INF)	return false;
		if((v	= props(4))==INF)	return false;
		rho	=	1/v;
		if((s	= props(5))==INF)	return false;
		if((u	= props(6))==INF)	return false;

		props = R22SCTR.getprop(P,T,props);
		if(!R22SCTR.goodprops())		return false;
		if((k	= props(3))==INF) return false;
		if((mu	= props(4))==INF) return false;
		if((Cp	= props(5))==INF) return false;

		break;
	case TWOPHASE:			//two-phase
		T	= Ts;
		v	= x*vv+(1-x)*vl;
		rho	= 1/v;
		s	= x*sv+(1-x)*sl;
		u	= x*uv+(1-x)*ul;
		k	= x*kv+(1-x)*kl;
		mu	= x*muv+(1-x)*mul;
		Cp	= x*Cpv+(1-x)*Cpl;
		break;
	case SUPERHEATED:			//superheated
		props = R22SHTH.getprop(P,h,props);
		if(!R22SHTH.goodprops())		return false;
		if((T	= props(3))==INF) return false;
		if((v	= props(4))==INF) return false;
		rho	= 1/v;
		if((s	= props(5))==INF) return false;
		if((u	= props(6))==INF) return false;
		props = R22SHTR.getprop(P,T,props);
		if(!R22SHTR.goodprops())		return false;
		if((k	= props(3))==INF) return false;
		if((mu	= props(4))==INF) return false;
		if((Cp	= props(5))==INF) return false;
		break;
	default:
		break;
	}//switch
	UPDATE = true;
	return true;
}

bool r22state::setstateP(double Press)
{
	Vector props(12);
	P=Press;

	//Determine all saturation properties at given pressure
	props = R222PTH.getprop(P,props);
	if(!R222PTH.goodprops())	return false;
	if((Ts = props(2))==INF) return false;
	if((hl = props(3))==INF) return false;
	if((hv = props(4))==INF) return false;
	if((vl = props(5))==INF) return false;
	if((vv = props(6))==INF) return false;
	rhol = 1/vl;
	rhov = 1/vv;
	if((sl = props(7))==INF) return false;
	if((sv = props(8))==INF) return false;
	if((ul = props(9))==INF) return false;
	if((uv = props(10))==INF) return false;
	props = R222PTR.getprop(P,props);
	if(!R222PTR.goodprops())	return false;
	if((kl = props(3))==INF) return false;
	if((kv = props(4))==INF) return false;
	if((mul= props(5))==INF) return false;
	if((muv= props(6))==INF) return false;
	if((Cpl= props(7))==INF) return false;
	if((Cpv= props(8))==INF) return false;

	return true;
}//r22state::setstateP

bool r22state::setstateT(double T)
{
	double Pguess = 500.0, Ttemp,dP = 100.0, Terr=1, Ttol=1e-2;
	Vector props(12);
	bool	errsign = true, olderrsign;

	do
	{
		props = R222PTH.getprop(Pguess,props);
		if(!R222PTH.goodprops())	return false;
		Ttemp = props(2);
		Terr = T - Ttemp;
		olderrsign = errsign;
		if(Terr>0)
		{
			errsign = true;
			Pguess += dP;
		}
		else
		{
			errsign = false;
			Pguess -= dP;
		}
		if(errsign!=olderrsign)
			dP /= 10.0;
	}while(fabs(Terr)>Ttol);

	P  = Pguess;
	props = R222PTH.getprop(P,props);
	if(!R222PTH.goodprops())	return false;
	if((Ts = props(2))==INF) return false;
	if((hl = props(3))==INF) return false;
	if((hv = props(4))==INF) return false;
	if((vl = props(5))==INF) return false;
	if((vv = props(6))==INF) return false;
	rhol = 1/vl;
	rhov = 1/vv;
	if((sl = props(7))==INF) return false;
	if((sv = props(8))==INF) return false;
	if((ul = props(9))==INF) return false;
	if((uv = props(10))==INF) return false;
	props = R222PTR.getprop(P,props);
	if(!R222PTR.goodprops())	return false;
	if((kl = props(3))==INF) return false;
	if((kv = props(4))==INF) return false;
	if((mul= props(5))==INF) return false;
	if((muv= props(6))==INF) return false;
	if((Cpl= props(7))==INF) return false;
	if((Cpv= props(8))==INF) return false;

	T = Ts;
	return true;
}//r22state::setstateT

bool	r22state::setstateHS(double h, double s)
{
	r22state	state;
	double	Pressure, Pg1, Pg2, sres1, sres2,temp, diff;
	Pg1 = 500.0;	Pg2 = 525.0;
	int iters=0;
	bool done = false;

	do
	{
		if(!state.setstate(Pg1,h))
		{
			errorlog.Add("r22state::setstateHS ","Could not set state with given h-s values.");
			return false;
		}
		sres1	=	s - state.gets();
		if(!state.setstate(Pg2,h))
		{
			errorlog.Add("r22state::setstateHS ","Could not set state with given h-s values.");
			return false;
		}
		sres2	=	s -	state.gets();

		diff = fabs((sres1-sres2)/s);
		if(diff>STOL)
		{
			temp	=	Pg1 - sres1*(Pg1 - Pg2)/(sres1 - sres2);
			Pg1		=	Pg2;
			Pg2		=	temp;
			iters++;
			if(iters>100)
			{
				errorlog.Add("r22state::setstateHS ","Convergence not achieved in 100 iterations.");
				return false;
			}
		}
		else
			done = true;
	}while(!done);
	Pressure = (Pg1+Pg2)/2;
	return this->setstate(Pressure,h);
}

bool	r22state::setstatePV(double P, double v)
{
	r22state state;

	double hmin, hmax;
	if(!state.setstateP(P))
	{
		errorlog.Add("r22state::setstatePV ","Could not determine two-phase properties with given pressure.");
		return false;
	}

	if(v<state.getvl())			{	hmax = state.gethl();		hmin = hmax - 20.0;		}
	else if(v>state.getvv())	{	hmin = state.gethv();		hmax = hmin + 50.0;		}
	else						{	hmin = state.gethl();		hmax = state.gethv();	}

	double residue, hmid, vmin, vmax;
	bool done = false;
	int iters = 0;
	while(!done)
	{
		//Search by bisection
		iters++;
		state.setstate(P,hmin);
		vmin = state.getv();
		state.setstate(P,hmax);
		vmax = state.getv();
		hmid = (hmin + hmax)/2.0;
		state.setstate(P,hmid);
		residue = v - state.getv();
		if(fabs(hmin-hmax)>1e-6)
		{
			if(residue>0)	
				hmin = hmid;
			else
				hmax = hmid;
		}
		else
			done = true;
	}//while !done

	return this->setstate(P,hmid);
}

bool	r22state::setstatePT(double P, double T)
{
	r22state state;
	double Tsat;
	if(!state.setstateP(P))
	{
		errorlog.Add("r22state::setstatePT ","Could not determine two-phase properties with given pressure.");
		return false;
	}
	//Check for two-phase condition
	Tsat = state.getTs();
	if((T-Tsat)<1e-6)
	{
		errorlog.Add("r22state::setstatePT ","Pressure and temperature are not independant.");
		return false;
	}
	//Not two-phase

	double hmin, hmax;
	if(T<Tsat)			{	hmax = state.gethl();	hmin = hmax - 20.0;	}
	else				{	hmin = state.gethv();	hmax = hmin + 50.0;	}

	double residue, hmid, Tmin, Tmax;
	bool done = false;
	int iters = 0;
	while(!done)
	{
		iters++;
		state.setstate(P,hmin);
		Tmin = state.getT();
		state.setstate(P,hmax);
		Tmax = state.getT();
		hmid = (hmin + hmax)/2.0;
		state.setstate(P,hmid);
		residue = T - state.getT();
		if(fabs(hmin-hmax)>1e-6)
		{
			if(residue>0)
				hmin = hmid;
			else
				hmax = hmid;
		}
		else
			done = true;
	}//while ! done
	return this->setstate(P,hmid);

}

/////////////////    R  500  //////////////////////////////////////////////////////////////////////

r500state::r500state()		{	UPDATE = false;	}

r500state::~r500state() {}

int		r500state::setphase(double enth, double hl, double hv)
{
	double	temp;
	int		temp2;

	temp2 = int(enth*1000.0);
	temp = (int(enth*1000.0))/1000.0;

	if(temp>hl&&enth<=hv)	//two phase
	{
		x = (enth-hl)/(hv-hl);
		return 2;
	}
	else if(temp<=hl)	//subcooled
	{
		x=-1.0;
		return 1;
	}
	else						//superheated
	{
		x=1.1;
		return 3;
	}
}//setphase and quality

bool r500state::setstateP(double Press)
{
	Vector props(10);
	P=Press;

	//Determine all saturation properties at given pressure
	props = BULBTPTH.getprop(P,props);
	if(!BULBTPTH.goodprops())	return false;
	Ts = props(2);
	hl = props(3);
	hv = props(4);
	vl = props(5);
	vv = props(6);
	rhol = 1/vl;
	rhov = 1/vv;
	sl = props(7);
	sv = props(8);
	ul = props(9);
	uv = props(10);
	T = Ts;

	return true;
}//r500state::setstateP

bool r500state::setstateT(double T)
{
	double Ptemp;
	Vector	props(10);
	double Pg1, Pg2, Terr1, Terr2;

	Pg1 = 500.0;
	Pg2 = 1500.0;

	do
	{
		props = BULBTPTH.getprop(Pg1,props);
		if(!BULBTPTH.goodprops())	return false;
		Terr1 = T - props(2);
		props = BULBTPTH.getprop(Pg2,props);
		if(!BULBTPTH.goodprops())	return false;
		Terr2 = T - props(2);

		Ptemp = Pg1 - (Pg1 - Pg2)*Terr1/(Terr1 - Terr2);
		Pg1 = Pg2;
		Pg2 = Ptemp;

	}while(fabs(Terr1-Terr2)>TTOL);

	P  = (Pg1+Pg2)/2;
	props = BULBTPTH.getprop(P,props);
	if(!BULBTPTH.goodprops())	return false;
	Ts = props(2);
	hl = props(3);
	hv = props(4);
	vl = props(5);
	vv = props(6);
	rhol = 1/vl;
	rhov = 1/vv;
	sl = props(7);
	sv = props(8);
	ul = props(9);
	uv = props(10);
	T = Ts;

	return true;
}//r500state::setstateT

bool loadtables(int refrigerant, bool BW)
{
	switch(refrigerant)
	{
	case R134A:
		if(!R134A2PTH.load2ptab("Properties\\R134a\\r134ahatpth.txt"))
		{
			errorlog.Add("loadtables() ","R134a two-phase thermodynamic properties could not be loaded");
			return false;
		}

		if(!R134A2PTR.load2ptab("Properties\\R134a\\r134ahatptr.txt"))
		{
			errorlog.Add("loadtables() ","R134a two-phase transport properties could not be loaded");	
			return false;
		}

		if(!R134ASCTH.load1ptab("Properties\\R134a\\r134ahascth.txt"))
		{
			errorlog.Add("loadtables() ","R134a subcooled thermodynamic properties could not be loaded");
			return false;
		}

		if(!R134ASCTR.load1ptab("Properties\\R134a\\r134ahasctr.txt"))
		{
			errorlog.Add("loadtables() ","R134a subcooled transport properties could not be loaded");
			return false;
		}

		if(!R134ASHTH.load1ptab("Properties\\R134a\\r134ahashth.txt"))
		{
			errorlog.Add("loadtables() ","R134a superheated thermodynamic properties could not be loaded");
			return false;
		}

		if(!R134ASHTR.load1ptab("Properties\\R134a\\r134ahashtr.txt"))
		{
			errorlog.Add("loadtables() ","R134a superheated transport properties could not be loaded");
			return false;
		}

		break;
	case R12:
		if(!R122PTH.load2ptab("Properties\\R12\\r12tpth.txt"))
		{
			errorlog.Add("loadtables() ","R12 two-phase thermodynamic properties could not be loaded");
			return false;
		}

		if(!R122PTR.load2ptab("Properties\\R12\\r12tptr.txt"))
		{
			errorlog.Add("loadtables() ","R12 two-phase transport properties could not be loaded");	
			return false;
		}

		if(!R12SCTH.load1ptab("Properties\\R12\\r12scth.txt"))
		{
			errorlog.Add("loadtables() ","R12 subcooled thermodynamic properties could not be loaded");
			return false;
		}

		if(!R12SCTR.load1ptab("Properties\\R12\\r12sctr.txt"))
		{
			errorlog.Add("loadtables() ","R12 subcooled transport properties could not be loaded");
			return false;
		}

		if(!R12SHTH.load1ptab("Properties\\R12\\r12shth.txt"))
		{
			errorlog.Add("loadtables() ","R12 superheated thermodynamic properties could not be loaded");
			return false;
		}

		if(!R12SHTR.load1ptab("Properties\\R12\\r12shtr.txt"))
		{
			errorlog.Add("loadtables() ","R12 superheated transport properties could not be loaded");
			return false;
		}

		break;
	case R22:
		if(!R222PTH.load2ptab("Properties\\R22\\r22tpth.txt"))
		{
			errorlog.Add("loadtables() ","R22 two-phase thermodynamic properties could not be loaded");
			return false;
		}

		if(!R222PTR.load2ptab("Properties\\R22\\r22tptr.txt"))
		{
			errorlog.Add("loadtables() ","R22 two-phase transport properties could not be loaded");	
			return false;
		}

		if(!R22SCTH.load1ptab("Properties\\R22\\r22scth.txt"))
		{
			errorlog.Add("loadtables() ","R22 subcooled thermodynamic properties could not be loaded");
			return false;
		}

		if(!R22SCTR.load1ptab("Properties\\R22\\r22sctr.txt"))
		{
			errorlog.Add("loadtables() ","R22 subcooled transport properties could not be loaded");
			return false;
		}

		if(!R22SHTH.load1ptab("Properties\\R22\\r22shth.txt"))
		{
			errorlog.Add("loadtables() ","R22 superheated thermodynamic properties could not be loaded");
			return false;
		}

		if(!R22SHTR.load1ptab("Properties\\R22\\r22shtr.txt"))
		{
			errorlog.Add("loadtables() ","R22 superheated transport properties could not be loaded");
			return false;
		}
		break;
	}

	if(BW)
	{
		if(!WATERLQ.load2ptab("Properties\\waterprops.txt"))
		{
			errorlog.Add("loadtables() ","Water transport properties could not be loaded");
			return false;
		}

		if(!BULBTPTH.load2ptab("Properties\\r500tpth.txt"))
		{
			errorlog.Add("R500 two-phase thermodynamic properties could not be loaded");
			return false;
		}
	}

	return true;
}//loadtables

