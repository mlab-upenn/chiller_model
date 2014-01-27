#include "StdAfx.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "Headers\Protos.h"
#include "Headers\Components.h"
#include "Headers\Matrix.h"
#include "Headers\Properties.h"
#include "Headers\NumRecipes.h"
#include "Headers\ErrorLog.h"

#define	SAMESIGN(X,Y)	((X*Y)>0.0)?true:false

extern ErrorLog	errorlog;
using namespace std;    

#define	MIN(X,Y)	(X<Y)?X:Y
#define MAX(X,Y)	(X>Y)?X:Y

bool	PROP_TABLES_LOADED = false;

extern two_phase_table	WATERLQ, BULBTPTH;

extern two_phase_table R134A2PTH, R134A2PTR;
extern one_phase_table R134ASCTH, R134ASCTR, R134ASHTH, R134ASHTR;

extern two_phase_table	R122PTH, R122PTR;
extern one_phase_table R12SCTH, R12SCTR, R12SHTH, R12SHTR;

extern two_phase_table	R222PTH, R222PTR;
extern one_phase_table R22SCTH, R22SCTR, R22SHTH, R22SHTR;

/******************************************************************************************/
/******************F I N I T E * V O L U M E * F O R M U L A T I O N***********************/
/****************************************O F***********************************************/
/***************S H E L L * - * T U B E * H E A T * E X C H A N G E R S********************/
/******************************************************************************************/

FVShellTubeHX::FVShellTubeHX()
{
	iGCOLS = 6;	iSCOLS = 5;
	VGeometry.make(10);
	VState.make(3);
	MState	=	new Vector [iSCOLS+1];

	iN_INPUTS	=	5;
	iN_OUTPUTS	=	9;
	VINPUTS.make(iN_INPUTS);
	VOUTPUTS.make(iN_OUTPUTS);

	dORIGREFQTY = dDYNREFQTY = dTOTMASSIN = dTOTMASSOUT = dREFMASSIMB = 0.0;

}

bool	FVShellTubeHX::fnDefineGeometry(const char* geomfile)
{
	ifstream	gfile(geomfile);
	if(!gfile)	return false;

	double tubeid, tubeod, tubelen, tubeCp, tuberho, foulingfactor, shellid;
	gfile >> iHXTYPE >> NNODES >> NTUBES;

	if(iHXTYPE!=EVAP&&iHXTYPE!=COND)
	{
		errorlog.Add("FVShellTubeHX::fnDefineGeometry ","Invalid heat exchanger type.");
		return false;
	}
	if(NNODES<=0||(int)NNODES!=NNODES)
	{
		errorlog.Add("FVShellTubeHX::fnDefineGeometry ","Invalid number of nodes.");
		return false;
	}
	if(NTUBES<=0||(int)NTUBES!=NTUBES)
	{
		errorlog.Add("FVShellTubeHX::fnDefineGeometry ","Invalid number of tubes.");
		return false;
	}

	gfile >> tubeid >> tubeod >> tubelen >> tubeCp >> tuberho >> foulingfactor >> shellid;
	if(tubeid>=tubeod||tubeod>=tubelen||tubeCp<=0||tuberho<=0||foulingfactor<0||foulingfactor>0.5||shellid<=0)
	{
		errorlog.Add("FVShellTubeHX::fnDefineGeometry ","Invalid tube dimensions and/or material properties.");
		return false;
	}

	VGeometry(1,tubeid);	VGeometry(2,tubeod);	VGeometry(3,tubelen);	VGeometry(4,foulingfactor);	VGeometry(5,shellid);

	double TotShellVolume, TotHTArea, TotTubeCap, TotWatCap, farea;
	double waterrho = 995.0, waterCp = 4.1868;

	farea			=	(MYPI/4)*(shellid*shellid-NTUBES*tubeod*tubeod);
	TotShellVolume	=	TotVolume = farea*tubelen;
	TotHTArea		=	NTUBES*MYPI*tubeod*tubelen;
	TotTubeCap		=	NTUBES*(MYPI/4)*(tubeod*tubeod-tubeid*tubeid)*tubelen*tuberho*tubeCp;
	TotWatCap		=	NTUBES*(MYPI/4)*tubeid*tubeid*tubelen*waterrho*waterCp;

	VGeometry(6,TotShellVolume/NNODES);
	VGeometry(7,farea);
	VGeometry(8,TotHTArea/NNODES);
	VGeometry(9,TotTubeCap/NNODES);
	VGeometry(10,TotWatCap/NNODES);

//enthalpies - tube temps - water temps - refrigerant flow rates - Nodal refrigerant side heat transfer rates
	for(int i=1;i<=iSCOLS;i++)
		MState[i].make(NNODES);

	gfile.close();

	A.make(NNODES+1,NNODES+1);
	B.make(NNODES+1);
	X.make(NNODES+1);
	ais.make(NNODES);
	bis.make(NNODES);
	cis.make(NNODES);
	dis.make(NNODES);
	partsumsa.make(NNODES+1);
	rpropertylist.make(20);
	wpropertylist.make(5);

	return true;

}

int		FVShellTubeHX::NODES()
{
	return	NNODES;
}

double	FVShellTubeHX::Volume()
{
	return	TotVolume;
}

double	FVShellTubeHX::fnWatHTCoeff(Vector& wpropertylist)
{

	double	Pr,Re,Nu,rho,mu,Cp,k,h,f,vel,tubeid;
	mu	=	wpropertylist(1);
	k	=	wpropertylist(2);
	Cp	=	wpropertylist(3);
	rho	=	wpropertylist(4);
	vel	=	wpropertylist(5);

	tubeid = VGeometry(1);

	Pr	= mu*Cp*1e3/k;
	Re	= rho*vel*tubeid/mu;

	switch(iHXTYPE)
	{
	case EVAP:
		f	=	-0.5087+0.2768*log10(Re)-0.0339*log10(Re)*log10(Re);
		Nu	=	((f/8)*(Re-1000)*Pr)/(1+12.7*sqrt(f/8)*(pow(Pr,0.67)-1));
		h	=	1.5*k*1e-3*Nu/tubeid;
		break;
	case COND:
		f	=	1/pow(0.79*log(Re)-1.64,2);
		Nu	=	(f*Re*Pr/8)/(1.07+12.7*sqrt(f/8)*(pow(Pr,0.67)-1));
		h	=	2.4*k*1e-3*Nu/tubeid;
		break;
	}//switch

	h	=	(iHXTYPE==COND)?3.0*h:1.25*h;
	return h;
}

double	FVShellTubeHX::fnRefHTCoeff(Vector& rpropertylist,double Tt)
{
		double Re, Pr, Nu, h, vel, tubeod,Qflux,qual, qualmin=0.05,qualmax=0.95;
		double Tsat, hfg, g=9.81,Cpl,rhol,rhov,kl,mul;
		double hfg1, den, num, Ja;
		double mur,Cpr,kr,rhor;
		double muv,Cpv,kv;
		int phase;
		double h2phase, hsuper;

		tubeod	= VGeometry(2);
		Tsat	= rpropertylist(1);
		mur		= rpropertylist(3);
		Cpr		= rpropertylist(4);
		kr		= rpropertylist(5);
		rhor	= rpropertylist(6);
		phase	= (int)rpropertylist(7);
		hfg		= rpropertylist(8)-rpropertylist(9);
		Cpl		= rpropertylist(10);
		rhol	= rpropertylist(11);
		rhov	= rpropertylist(12);
		kl		= rpropertylist(13);
		mul		= rpropertylist(14);
		vel		= rpropertylist(15);
		Qflux	= rpropertylist(16);
		qual	= rpropertylist(17);
		muv		= rpropertylist(18);
		Cpv		= rpropertylist(19);
		kv		= rpropertylist(20);

		switch(phase)
		{
		case SUBCOOLED:
			Pr	=	(mur*Cpr*1e3)/kr;
			Re	=	(rhor*vel*tubeod)/mur;
			Nu	=	1.13*0.5*pow(Re,0.5)*Pr;
			h	=	(Nu*kr/tubeod)*1e-3;		//kW/m^2-K
			break;
		case TWOPHASE:
			switch(iHXTYPE)
			{
			case COND:
				double	hsubcooled;
				Ja		=	Cpl*fabs(Tsat-Tt)/hfg;
				hfg1	=	hfg*(1+0.68*Ja)*1e3;
				num		=	g*rhol*(rhol-rhov)*pow(kl,3)*hfg1;
				den		=	mul*fabs(Tsat-Tt)*tubeod;	
				if(den>0.0&&num>0.0)
					h2phase	=	7.5*0.729*pow(num/den,0.25)*1e-3;	//kW/m^2-K
				else
					h2phase = 0.0;
				if(qual<qualmin)
				{
					double	Pr0, Re0, Nu0;
					Pr0	=	(mul*Cpl*1e3)/kl;
					Re0	=	(rhol*vel*tubeod)/mul;
					Nu0	=	1.13*0.5*pow(Re0,0.5)*Pr0;
					hsubcooled	=	(Nu0*kl/tubeod)*1e-3;
					h	=	hsubcooled + (qual/qualmin)*(h2phase-hsubcooled);
				}
				else if(qual > qualmax)
				{
					double Pr1, Re1, Nu1;
					Pr1	=	(muv*Cpv*1e3)/kr;
					Re1	=	(rhov*vel*tubeod)/muv;
					Nu1	=	1.13*pow(Re1,0.5)*Pr1;
					hsuper=	2*(Nu1*kv/tubeod)*1e-3;
					h	=	hsuper + ((1-qual)/(1-qualmax))*(h2phase-hsuper);
				}
				else
					h	=	h2phase;
				break;
			case EVAP:
				h2phase =	12.492 + Qflux/2.22;	//kW/m^2-K
				if(qual<=qualmax)
					h	=	h2phase;
				else
				{
					double Pr0, Re0, Nu0;
					Pr0	= (muv*Cpv*1e3)/kv;
					Re0	= (rhov*vel*tubeod)/muv;
					Nu0 = 1.13*0.5*pow(Re0,0.5)*Pr0;
					hsuper = 4*(Nu0*kv/tubeod)*1e-3;
					h	=	hsuper + ((1-qual)/(1-qualmax))*(h2phase-hsuper);
				}
				break;
			}//switch iHXTYPE
			break;
		case SUPERHEATED:
			Pr	=	(mur*Cpr*1e3)/kr;
			Re	=	(rhor*vel*tubeod)/mur;
			Nu	=	1.13*pow(Re,0.5)*Pr;
			h	=	2*(Nu*kr/tubeod)*1e-3;
			break;
		}//switch phase

		h	=	(iHXTYPE==COND)?3.0*h:1.25*h;
		return h;
}

void	FVShellTubeHX::fnGetState(double& P, Matrix& profiles)
{
	P = VState(1);

	int nnodes;
	nnodes = NNODES;
	profiles.make(nnodes,iSCOLS);

	for(int i=1;i<=iSCOLS;i++)
		profiles.putcol(i,MState[i]);

}

void	FVShellTubeHX::fnGetState(Vector& sv1)
{
	sv1 = VState;
}

bool	FVShellTubeHX::fnSetState(double time, double svs1,Matrix& svs2)
{

	dTIME	=	time;
	char errstr[512];

	int nnodes;
	nnodes = NNODES;
	if(svs2.rows()!=nnodes||svs2.cols()!=iSCOLS)
	{
		errorlog.Add("FVShellTubeHX::fnSetState ","State and Initialization matrices size mismatch.");
		return false;
	}

	for(int i=1;i<=iSCOLS;i++)
		svs2.getcol(i,MState[i]);

	r134astate	exitstate;

	VState(1,svs1);					//Pressure
	VState(2,MState[1](nnodes));	//Exit enthalpy
	VState(3,MState[3](1));			//Leaving water temperature

	double P,RefMass=0.0, vol;
	Vector h;
	h.make(nnodes);
	P	=	svs1;
	h	=	MState[1];
	vol =	VGeometry(6);
	r134astate	state;
	for(int i=1;i<=nnodes;i++)
	{
		if(!state.setstate(P,h(i)))	
		{
			sprintf(errstr,"Could not set refrigerant properties with %fkPa and %fkJ/kg",P,h(i));
			errorlog.Add("FVShellTubeHX::fnSetState ",errstr);
			return false;
		}
		RefMass	+=	vol*state.getrho();
	}

	dORIGREFQTY	= RefMass;
	dDYNREFQTY	= RefMass;

	return true;
}

bool	FVShellTubeHX::fnSetState(double time, double Pressure, double RefMass)
{
	dTIME	=	time;
	char errstr[512];

	int nnodes;
	nnodes = NNODES;

	dORIGREFQTY	= RefMass;
	dDYNREFQTY	= RefMass;

	//Compute refrigerant state from given pressure and total mass
	double v;
	v = TotVolume/RefMass;
	r134astate state;
	if(!state.setstatePV(Pressure,v))
	{
		sprintf(errstr,"Could not set refrigerant properties with %fkPa and %fm^3/kg",Pressure,v);
		errorlog.Add("FVShellTubeHX::fnSetState ",errstr);
		return false;
	}

	//Set enthalpy distribution and tube/water temperatures
	double h, T;
	h = state.geth();
	T = state.getT();
	for(int j=1;j<=nnodes;j++)
	{
		MState[1](j,h);
		MState[2](j,T);
		MState[3](j,T);
		MState[4](j,0.0);
		MState[5](j,0.0);
	}

	//Set leaving water temperature
	VState(1,Pressure);		//Pressure
	VState(2,h);			//Exit enthalpy
	VState(3,T);			//Leaving water temperature

	return true;
}

bool	FVShellTubeHX::fnPopulateAB(double mrin, double hrin, double mrout,Vector& hi, Vector& Qri)
{
	int nnodes;
	nnodes = NNODES;
	double temp;
	for(int i=1;i<=nnodes+1;i++)
	{
		if(i==1)
			B(i,mrin-mrout);
		else if(i==2)
			B(i,mrin*hrin - mrout*hi(nnodes) - sumQ);
		else
		{
			temp = (i==3)?hrin:hi(i-3);
			B(i,mrin*(temp - hi(i-2))-Qri(i-2));
		}
		for(int j=1;j<=nnodes+1;j++)
		{
			if(i==1)
			{
				temp = (j==1)?partsumsa(nnodes+1):bis(j-1);
				A(i,j,temp);
			}
			else if(i==2)
			{
				temp = (j==1)?sumc:dis(j-1);
				A(i,j,temp);
			}
			else
			{
				if(j==1)
				{
					temp = (i==3)?hrin:hi(i-3);
					A(i,j,cis(i-2)-ais(i-2)*hi(i-2) - partsumsa(i-2)*(hi(i-2)-temp));
				}
				else if (j==i-1)
					A(i,j,dis(i-2)-bis(i-2)*hi(i-2));
				else if (j>1 && j<i-1)
				{
					temp = (i==3)?hrin:hi(i-3);
					A(i,j,-bis(j-1)*(hi(i-2)-temp));
				}
				else
					A(i,j,0.0);
			}
		}
	}
	return true;
}

bool	FVShellTubeHX::fnStateDerivs(double P,Vector& h,Vector& Tt, Vector& Tw, Vector& mr0, Vector& Qr0,
									 double& Pdot, Vector& hdot, Vector& Ttdot, Vector& Twdot,
									 double& M, double& vrhou)
{
	//Capture geometry
	int nnodes;
	nnodes = NNODES;
	double	tubeid, tubeod, diaratio, foulfactor;
	tubeid	=	VGeometry(1);	tubeod	=	VGeometry(2);	diaratio = tubeid/tubeod;
	foulfactor = 1.0 - VGeometry(4);
	double	vol, farea, htarea, tcap, wcap;
	vol		=	VGeometry(6);
	farea	=	VGeometry(7);
	htarea	=	VGeometry(8);
	tcap	=	VGeometry(9);
	wcap	=	VGeometry(10);

	r134astate	state;
	double mi_1,mi;
	double Twi_1;

	//Capture boundary conditions
	double mwat, Twin, mrin, mrout,hrin;
	mwat = VINPUTS(1);	Twin = VINPUTS(2);	mrin = VINPUTS(3);	mrout = VINPUTS(4);	hrin = VINPUTS(5);

	//Scratch variables
	watstate water;
	double	whtcoeff, rhtcoeff, Tr, Qw, Cpw;
	double	a,b,c,d;
	Vector	x;		x.make(nnodes);
	Vector Qflux;	Qflux.make(nnodes);
	double	rho;

	bool	done = false;
	int		iters = 0;

	double area, vel;
	area	= MYPI*tubeid*tubeid/4.0;
	water.setstate(Twin);
	vel		=	2*(mwat/NTUBES)/(water.getrho()*area);
	wpropertylist(1,water.getmu());
	wpropertylist(2,water.getk());
	wpropertylist(3,water.getCp());
	wpropertylist(4,water.getrho());
	wpropertylist(5,vel);
	whtcoeff = foulfactor*fnWatHTCoeff(wpropertylist);
	Cpw	=	water.getCp();

	vrhou = 0.0;
	M = 0.0;

	partsumsa.reset();
	sumc = 0.0;
	sumQ = 0.0;
	for(int i=1;i<=nnodes;i++)
	{
		state.setstate(P,h(i));
		x(i,state.getx());
		rho	=	state.getrho();
		Tr	=	state.getT();
		M	+=	vol*rho;
		vrhou	+=	vol*rho*state.getu();
		mi_1	=	(i==i)?mrin:mr0(i-1);
		mi		=	(i==nnodes)?mrout:mr0(i);
		vel		=	fabs(mi_1+mi)/(2*farea*rho);
		Qflux(i,fabs(Qr0(i)/htarea));
		rpropertylist(1,state.getTs());		rpropertylist(3,state.getmu());		rpropertylist(4,state.getCp());
		rpropertylist(5,state.getk());		rpropertylist(6,state.getrho());	rpropertylist(7,state.getphase());
		rpropertylist(8,state.gethv());		rpropertylist(9,state.gethl());		rpropertylist(10,state.getCpl());
		rpropertylist(11,state.getrhol());	rpropertylist(12,state.getrhov());	rpropertylist(13,state.getkl());
		rpropertylist(14,state.getmul());	rpropertylist(15,vel);				rpropertylist(16,Qflux(i));
		rpropertylist(17,x(i));				rpropertylist(18,state.getmuv());	rpropertylist(19,state.getCpv());
		rpropertylist(20,state.getkv());

		rhtcoeff	=	fnRefHTCoeff(rpropertylist,Tt(i));
		Qr0(i,rhtcoeff*htarea*(Tr-Tt(i)));
		Qw	=	whtcoeff*htarea*diaratio*(Tt(i)-Tw(i));
		Ttdot(i,(Qr0(i)-Qw)/tcap);
		Twi_1	=	(i==nnodes)?Twin:Tw(i+1);
		Twdot(i,(mwat*Cpw*(Twi_1-Tw(i))+Qw)/wcap);
		a	=	vol*state.getdrdP();
		b	=	vol*state.getdrdh();
		c	=	vol*(h(i)*state.getdrdP()-1.0);
		d	=	vol*(h(i)*state.getdrdh()+state.getrho());
		ais(i,a);	bis(i,b);	cis(i,c);	dis(i,d);
		partsumsa(i+1,partsumsa(i)+a);
		sumc += c;
		sumQ += Qr0(i);
	}//for i
	fnPopulateAB(mrin,hrin,mrout,h,Qr0);
	solveLU(A,B,X);
	Pdot = X(1);
	for(int i=1;i<=nnodes;i++)
	{
		hdot(i,X(i+1));
		if(i==1)
			mr0(i,mrin-ais(i)*Pdot-bis(i)*hdot(i));
		else if(i==nnodes)
			mr0(i,mrout);
		else
			mr0(i,mr0(i-1)-ais(i)*Pdot-bis(i)*hdot(i));
	}

	return true;
}

bool	FVShellTubeHX::fnReEvalStateExpE(double& newstate1, Matrix& newstate2, double& HXEnBal)
{
	//Capture geometry
	int nnodes;
	nnodes = NNODES;
	double vol;
	vol = VGeometry(6);

	//Capture state
	double P;
	P	=	VState(1);
	Vector hi;		hi.make(nnodes);	hi	=	MState[1];
	Vector Tt;		Tt.make(nnodes);	Tt	=	MState[2];
	Vector Tw;		Tw.make(nnodes);	Tw	=	MState[3];
	Vector mr0;		mr0.make(nnodes);	mr0	=	MState[4];
	Vector Qr0;		Qr0.make(nnodes);	Qr0 =	MState[5];

	//Capture boundary conditions
	double mwat, Twin, mrin, mrout,hrin;
	mwat = VINPUTS(1);	Twin = VINPUTS(2);	mrin = VINPUTS(3);	mrout = VINPUTS(4);	hrin = VINPUTS(5);

	//Integration variables
	double time=0, tstart=0.0, tend=1.0, tstep;
	double dPdt;
	Vector dhdt, dTtdt, dTwdt, dTshdt;
	dhdt.make(nnodes);	dTtdt.make(nnodes);	dTwdt.make(nnodes);	dTshdt.make(nnodes);

	//Scratch variables
	double dP1;
	Vector dh1, dh2;	dh1.make(nnodes);
	Vector dTt1,dTt2;	dTt1.make(nnodes);
	Vector dTw1,dTw2;	dTw1.make(nnodes);

	double Porig;
	Vector horig;		horig.make(nnodes);
	Vector Ttorig;		Ttorig.make(nnodes);
	Vector Tworig;		Tworig.make(nnodes);
	double	Cpw;

	r134astate	state;
	double M1, M2;

	Porig = P;
	horig = hi;
	Ttorig = Tt;
	Tworig = Tw;
	tstart = 0.0;

	double VRHOU_0=0.0, VRHOU_1=0.0;
	bool	done = false;
	int		iters = 0;
	char errstr[512];
	while(tstart<tend)
	{
		tstep	=	((tstart+dTSTEP)>tend)?(tend-tstart):dTSTEP;
		fnStateDerivs(P,hi,Tt,Tw,mr0,Qr0,dPdt,dhdt,dTtdt,dTwdt,M1,VRHOU_1);
		if(tstart==0.0)	VRHOU_0 = VRHOU_1;
		dP1		=	dPdt*tstep;
		P		=	Porig +	dP1;
		for(int i=1;i<=nnodes;i++)
		{
			dh1(i,dhdt(i)*tstep);
			dTt1(i,dTtdt(i)*tstep);
			dTw1(i,dTwdt(i)*tstep);
			hi(i,horig(i)+dh1(i));
			Tt(i,Ttorig(i)+dTt1(i));
			Tw(i,Tworig(i)+dTw1(i));
		}//End of Integration
		M2	=	0.0;
		VRHOU_1	= 0.0;
		for(int i=1;i<=nnodes;i++)
		{
			if(!state.setstate(P,hi(i)))
			{
				sprintf(errstr,"Cannot set properties at %fkPa and %fkJ/kg",P,hi(i));
				errorlog.Add("FVShellTubeHX::fnReEvalState ",errstr);
			}
			M2	+=	vol*state.getrho();
			VRHOU_1	+=	vol*state.getrho()*state.getu();
		}
		dDYNREFQTY	=	M2;
		dTOTMASSIN	+=	mrin*tstep;
		dTOTMASSOUT	+=	mrout*tstep;
		dREFMASSIMB	=	(dDYNREFQTY - dORIGREFQTY)-(dTOTMASSIN-dTOTMASSOUT);

		tstart	+=	tstep;
		Porig = P;
		horig = hi;
		Ttorig= Tt;
		Tworig= Tw;
	}//while tstart < tend

	//Energy Balance computation for the 1 sec step;
	double Pinit;
	Vector hinit;		hinit.make(nnodes);
	Vector Ttinit;		Ttinit.make(nnodes);
	Vector Twinit;		Twinit.make(nnodes);
	Pinit	=	VState(1);
	hinit	=	MState[1];
	Ttinit	=	MState[2];
	Twinit	=	MState[3];

	double tcap, wcap;
	tcap = VGeometry(9);
	wcap = VGeometry(10);

	double	Qrin, Qrout, Qwin, Qwout, Qrstored, Qtstored, Qwstored, Qtotstored;
	Cpw		=	4.1868;	//kJ/kg-K
	Qrin	=	mrin*hrin;		Qrout	=	mrout*hinit(nnodes);
	Qwin	=	mwat*Cpw*Twin;	Qwout	=	mwat*Cpw*Twinit(1);
	Qrstored=	VRHOU_1-VRHOU_0;
	Qtstored = 0.0;
	Qwstored = 0.0;
	for(int i=1;i<=nnodes;i++)
	{
		Qtstored += tcap*(Tt(i)-Ttinit(i));
		Qwstored += wcap*(Tw(i)-Twinit(i));
	}
	Qtotstored	=	Qrstored + Qtstored + Qwstored;
	HXEnBal	=	(Qrin + Qwin) - (Qrout+Qwout) - Qtotstored;
	//End of Energy Balance computation

	newstate1 = P;
	newstate2.putcol(1,hi);
	newstate2.putcol(2,Tt);
	newstate2.putcol(3,Tw);
	newstate2.putcol(4,mr0);
	newstate2.putcol(5,Qr0);

	return true;
}

bool	FVShellTubeHX::fnReEvalStateExpEPC(double& newstate1, Matrix& newstate2, double& HXEnBal)
{
	//Capture geometry
	int nnodes;
	nnodes = NNODES;
	double vol;
	vol = VGeometry(6);

	//Capture state
	double P;
	P	=	VState(1);
	Vector hi;		hi.make(nnodes);	hi	=	MState[1];
	Vector Tt;		Tt.make(nnodes);	Tt	=	MState[2];
	Vector Tw;		Tw.make(nnodes);	Tw	=	MState[3];
	Vector mr0;		mr0.make(nnodes);	mr0	=	MState[4];
	Vector Qr0;		Qr0.make(nnodes);	Qr0 =	MState[5];

	//Capture boundary conditions
	double mwat, Twin, mrin, mrout,hrin;
	mwat = VINPUTS(1);	Twin = VINPUTS(2);	mrin = VINPUTS(3);	mrout = VINPUTS(4);	hrin = VINPUTS(5);

	//Integration variables
	double time=0, tstart=0.0, tend=1.0, tstep;
	double dPdt;
	Vector dhdt, dTtdt, dTwdt, dTshdt;
	dhdt.make(nnodes);	dTtdt.make(nnodes);	dTwdt.make(nnodes);	dTshdt.make(nnodes);

	//Scratch variables
	double dP1, dP2;
	Vector dh1, dh2;	dh1.make(nnodes);	dh2.make(nnodes);
	Vector dTt1,dTt2;	dTt1.make(nnodes);	dTt2.make(nnodes);
	Vector dTw1,dTw2;	dTw1.make(nnodes);	dTw2.make(nnodes);

	double Porig;
	Vector horig;		horig.make(nnodes);
	Vector Ttorig;		Ttorig.make(nnodes);
	Vector Tworig;		Tworig.make(nnodes);
	double	Cpw;

	r134astate	state;
	double M1, M2;

	Porig = P;
	horig = hi;
	Ttorig = Tt;
	Tworig = Tw;
	tstart = 0.0;

	double VRHOU_0=0.0, VRHOU_1=0.0;
	bool	done = false;
	int		iters = 0;
	char errstr[512];
	while(tstart<tend)
	{
		tstep	=	((tstart+dTSTEP)>tend)?(tend-tstart):dTSTEP;
		{
			//Predictor step
			fnStateDerivs(P,hi,Tt,Tw,mr0,Qr0,dPdt,dhdt,dTtdt,dTwdt,M1,VRHOU_1);
			if(tstart==0.0)	VRHOU_0 = VRHOU_1;
		}
		done = false;
		while(!done)
		{
			dP1		=	dPdt*tstep;
			P		=	Porig +	dP1;
			for(int i=1;i<=nnodes;i++)
			{
				dh1(i,dhdt(i)*tstep);
				dTt1(i,dTtdt(i)*tstep);
				dTw1(i,dTwdt(i)*tstep);
				hi(i,horig(i)+dh1(i));
				Tt(i,Ttorig(i)+dTt1(i));
				Tw(i,Tworig(i)+dTw1(i));

			}//End of Predictor step
			{
				//Corrector step
				fnStateDerivs(P,hi,Tt,Tw,mr0,Qr0,dPdt,dhdt,dTtdt,dTwdt,M2,VRHOU_1);
			}
			dP2		=	dPdt*tstep;
			for(int i=1;i<=nnodes;i++)
			{
				dh2(i,dhdt(i)*tstep);
				dTt2(i,dTtdt(i)*tstep);
				dTw2(i,dTwdt(i)*tstep);
			}

			P	=	Porig + (dP1+dP2)/2.0;

			M2	=	0.0;
			VRHOU_1	= 0.0;
			for(int i=1;i<=nnodes;i++)
			{
				hi(i,horig(i)+(dh1(i)+dh2(i))/2.0);
				Tt(i,Ttorig(i)+(dTt1(i)+dTt2(i))/2.0);
				Tw(i,Tworig(i)+(dTw1(i)+dTw2(i))/2.0);

				if(!state.setstate(P,hi(i)))
				{
					sprintf(errstr,"Cannot set properties at %fkPa and %fkJ/kg",P,hi(i));
					errorlog.Add("FVShellTubeHX::fnReEvalState ",errstr);
				}
				M2	+=	vol*state.getrho();
				VRHOU_1	+=	vol*state.getrho()*state.getu();
			}
			if(bTSIZING)
			{
				double IMB;
				IMB = fabs((M2 - dORIGREFQTY) - (dTOTMASSIN+mrin*tstep - dTOTMASSOUT-mrout*tstep));
				if(IMB > MAXIMB && iters < MAXITERS)
				{
					iters++;
					tstep	= tstep/2.0;
					tstep	= (tstep<MINDT)?MINDT:tstep;
				}
				else if(IMB < MINIMB)
				{
					tstep	+=	DTINC;
					tstep = (tstep>MAXDT)?MAXDT:tstep;
					done = true;
				}
				else
					done = true;
			}
			else
				done = true;

			done = true;
			if(done)
			{
				iters	=	0;
				dTSTEP	=	tstep;
				dDYNREFQTY	=	M2;
				dTOTMASSIN	+=	mrin*tstep;
				dTOTMASSOUT	+=	mrout*tstep;
				dREFMASSIMB	=	(dDYNREFQTY - dORIGREFQTY)-(dTOTMASSIN-dTOTMASSOUT);
			}
		}//while(!done)...End of Predictor-Corrector

		tstart	+=	tstep;
		Porig = P;
		horig = hi;
		Ttorig= Tt;
		Tworig= Tw;
	}//while tstart < tend

	//Energy Balance computation for the 1 sec step;
	double Pinit;
	Vector hinit;		hinit.make(nnodes);
	Vector Ttinit;		Ttinit.make(nnodes);
	Vector Twinit;		Twinit.make(nnodes);
	Pinit	=	VState(1);
	hinit	=	MState[1];
	Ttinit	=	MState[2];
	Twinit	=	MState[3];

	double tcap, wcap;
	tcap = VGeometry(9);
	wcap = VGeometry(10);

	double	Qrin, Qrout, Qwin, Qwout, Qrstored, Qtstored, Qwstored, Qtotstored;
	Cpw		=	4.1868;	//kJ/kg-K
	Qrin	=	mrin*hrin;		Qrout	=	mrout*hinit(nnodes);
	Qwin	=	mwat*Cpw*Twin;	Qwout	=	mwat*Cpw*Twinit(1);
	Qrstored=	VRHOU_1-VRHOU_0;
	Qtstored = 0.0;
	Qwstored = 0.0;
	for(int i=1;i<=nnodes;i++)
	{
		Qtstored += tcap*(Tt(i)-Ttinit(i));
		Qwstored += wcap*(Tw(i)-Twinit(i));
	}
	Qtotstored	=	Qrstored + Qtstored + Qwstored;
	HXEnBal	=	(Qrin + Qwin) - (Qrout+Qwout) - Qtotstored;
	//End of Energy Balance computation

	newstate1 = P;
	newstate2.putcol(1,hi);
	newstate2.putcol(2,Tt);
	newstate2.putcol(3,Tw);
	newstate2.putcol(4,mr0);
	newstate2.putcol(5,Qr0);

	return true;
}

bool	FVShellTubeHX::fnReEvalStateImpEPC(double& newstate1, Matrix& newstate2, double& HXEnBal)
{
	//Capture geometry
	int nnodes;
	nnodes = NNODES;
	double vol;
	vol = VGeometry(6);


	//Capture state
	double P;
	P	=	VState(1);
	Vector hi;		hi.make(nnodes);	hi	=	MState[1];
	Vector Tt;		Tt.make(nnodes);	Tt	=	MState[2];
	Vector Tw;		Tw.make(nnodes);	Tw	=	MState[3];
	Vector mr0;		mr0.make(nnodes);	mr0	=	MState[4];
	Vector Qr0;		Qr0.make(nnodes);	Qr0 =	MState[5];

	//Capture boundary conditions
	double mwat, Twin, mrin, mrout,hrin;
	mwat = VINPUTS(1);	Twin = VINPUTS(2);	mrin = VINPUTS(3);	mrout = VINPUTS(4);	hrin = VINPUTS(5);

	//Integration variables
	double time=0, tstart=0.0, tend=1.0, tstep;
	double dPdt;
	Vector dhdt, dTtdt, dTwdt, dTshdt;
	dhdt.make(nnodes);	dTtdt.make(nnodes);	dTwdt.make(nnodes);	dTshdt.make(nnodes);

	//Scratch variables
	double	Pres;
	Vector	hres;	hres.make(nnodes);
	Vector	Ttres;	Ttres.make(nnodes);
	Vector	Twres;	Twres.make(nnodes);
	Vector	Residuals;	Residuals.make(3*nnodes+1);

	double Porig;
	Vector horig;		horig.make(nnodes);
	Vector Ttorig;		Ttorig.make(nnodes);
	Vector Tworig;		Tworig.make(nnodes);
	double	Cpw;

	r134astate	state;
	double M2;

	Porig = P;
	horig = hi;
	Ttorig = Tt;
	Tworig = Tw;
	tstart = 0.0;

	double VRHOU_0=0.0, VRHOU_1=0.0;
	bool	done = false;
	int		iters = 0;
	double	URF = 0.04;
//	char errstr[512];

	while(tstart<tend)
	{
//		tstep	=	((tstart+dTSTEP)>tend)?(tend-tstart):dTSTEP;
		tstep	=	1.0;
		bool done = false;
		while(!done)
		{
			fnStateDerivs(P,hi,Tt,Tw,mr0,Qr0,dPdt,dhdt,dTtdt,dTwdt,M2,VRHOU_1);
			Pres = Porig - (P - dPdt*tstep);
			Residuals(1,Pres);
			P	+= URF*Pres;
			for(int i=1;i<=nnodes;i++)
			{
				hres(i,horig(i) - (hi(i) - dhdt(i)*tstep));
				Residuals(3*i-1,hres(i));
				Ttres(i,Ttorig(i) - (Tt(i) - dTtdt(i)*tstep));
				Residuals(3*i,Ttres(i));
				Twres(i,Tworig(i) - (Tw(i) - dTwdt(i)*tstep));
				Residuals(3*i+1,Twres(i));

				hi(i,hi(i) + URF*hres(i));
				Tt(i,Tt(i) + URF*Ttres(i));
				Tw(i,Tw(i) + URF*Twres(i));
			}
			double test;
			test = Residuals.fmax();
			if(test < 1e-2)
				done = true;
		}//while(!done)...End of Predictor-Corrector

		tstart	+=	tstep;
		Porig = P;
		horig = hi;
		Ttorig= Tt;
		Tworig= Tw;
	}//while tstart < tend

	//Energy Balance computation for the 1 sec step;
	double Pinit;
	Vector hinit;		hinit.make(nnodes);
	Vector Ttinit;		Ttinit.make(nnodes);
	Vector Twinit;		Twinit.make(nnodes);
	Pinit	=	VState(1);
	hinit	=	MState[1];
	Ttinit	=	MState[2];
	Twinit	=	MState[3];

	double tcap, wcap;
	tcap = VGeometry(9);
	wcap = VGeometry(10);

	double	Qrin, Qrout, Qwin, Qwout, Qrstored, Qtstored, Qwstored, Qtotstored;
	Cpw		=	4.1868;	//kJ/kg-K
	Qrin	=	mrin*hrin;		Qrout	=	mrout*hinit(nnodes);
	Qwin	=	mwat*Cpw*Twin;	Qwout	=	mwat*Cpw*Twinit(1);
	Qrstored=	VRHOU_1-VRHOU_0;
	Qtstored = 0.0;
	Qwstored = 0.0;
	for(int i=1;i<=nnodes;i++)
	{
		Qtstored += tcap*(Tt(i)-Ttinit(i));
		Qwstored += wcap*(Tw(i)-Twinit(i));
	}
	Qtotstored	=	Qrstored + Qtstored + Qwstored;
	HXEnBal	=	(Qrin + Qwin) - (Qrout+Qwout) - Qtotstored;
	//End of Energy Balance computation

	newstate1 = P;
	newstate2.putcol(1,hi);
	newstate2.putcol(2,Tt);
	newstate2.putcol(3,Tw);
	newstate2.putcol(4,mr0);
	newstate2.putcol(5,Qr0);

	return true;
}

bool	FVShellTubeHX::fnReEvalStateRK4(double& newstate1, Matrix& newstate2, double& HXEnBal)
{
	//Capture geometry
	int nnodes;
	nnodes = NNODES;
	double vol;
	vol = VGeometry(6);

	//Capture state
	double P;
	P	=	VState(1);
	Vector hi;		hi.make(nnodes);	hi	=	MState[1];
	Vector Tt;		Tt.make(nnodes);	Tt	=	MState[2];
	Vector Tw;		Tw.make(nnodes);	Tw	=	MState[3];
	Vector mr0;		mr0.make(nnodes);	mr0	=	MState[4];
	Vector Qr0;		Qr0.make(nnodes);	Qr0 =	MState[5];

	//Capture boundary conditions
	double mwat, Twin, mrin, mrout,hrin;
	mwat = VINPUTS(1);	Twin = VINPUTS(2);	mrin = VINPUTS(3);	mrout = VINPUTS(4);	hrin = VINPUTS(5);

	//Integration variables
	double time=0, tstart=0.0, tend=1.0, tstep;
	double dPdt;
	Vector dhdt, dTtdt, dTwdt, dTshdt;
	dhdt.make(nnodes);	dTtdt.make(nnodes);	dTwdt.make(nnodes);	dTshdt.make(nnodes);

	//Scratch variables
	double dP1, dP2, dP3, dP4;
	Vector dh1, dh2, dh3, dh4;
	dh1.make(nnodes);	dh2.make(nnodes);	dh3.make(nnodes);	dh4.make(nnodes);
	Vector dTt1, dTt2, dTt3, dTt4;
	dTt1.make(nnodes);	dTt2.make(nnodes);	dTt3.make(nnodes);	dTt4.make(nnodes);
	Vector dTw1, dTw2, dTw3, dTw4;
	dTw1.make(nnodes);	dTw2.make(nnodes);	dTw3.make(nnodes);	dTw4.make(nnodes);

	double Porig;	Porig = P;
	Vector horig;		horig.make(nnodes);		horig = hi;
	Vector Ttorig;		Ttorig.make(nnodes);	Ttorig = Tt;
	Vector Tworig;		Tworig.make(nnodes);	Tworig = Tw;
	double	Cpw;

	r134astate	state;
	double M1, M2;
	tstart = 0.0;

	double VRHOU_0=0.0, VRHOU_1=0.0;
	bool	done = false;
	int		iters = 0;
	char errstr[512];
	
	while(tstart<tend)
	{
		tstep	=	((tstart+dTSTEP)>tend)?(tend-tstart):dTSTEP;
		fnStateDerivs(P,hi,Tt,Tw,mr0,Qr0,dPdt,dhdt,dTtdt,dTwdt,M1,VRHOU_1);
		if(tstart==0.0)	VRHOU_0 = VRHOU_1;
		dP1		=	dPdt*tstep;
		P		=	Porig +	dP1/2.0;
		for(int i=1;i<=nnodes;i++)
		{
			dh1(i,dhdt(i)*tstep);
			dTt1(i,dTtdt(i)*tstep);
			dTw1(i,dTwdt(i)*tstep);
			hi(i,horig(i)+dh1(i)/2.0);
			Tt(i,Ttorig(i)+dTt1(i)/2.0);
			Tw(i,Tworig(i)+dTw1(i)/2.0);
		}
		fnStateDerivs(P,hi,Tt,Tw,mr0,Qr0,dPdt,dhdt,dTtdt,dTwdt,M1,VRHOU_1);
		dP2		=	dPdt*tstep;
		P		=	Porig + dP2/2.0;
		for(int i=1;i<=nnodes;i++)
		{
			dh2(i,dhdt(i)*tstep);
			dTt2(i,dTtdt(i)*tstep);
			dTw2(i,dTwdt(i)*tstep);
			hi(i,horig(i)+dh2(i)/2.0);
			Tt(i,Ttorig(i)+dTt2(i)/2.0);
			Tw(i,Tworig(i)+dTw2(i)/2.0);
		}
		fnStateDerivs(P,hi,Tt,Tw,mr0,Qr0,dPdt,dhdt,dTtdt,dTwdt,M1,VRHOU_1);
		dP3		=	dPdt*tstep;
		P		=	Porig + dP3;
		for(int i=1;i<=nnodes;i++)
		{
			dh3(i,dhdt(i)*tstep);
			dTt3(i,dTtdt(i)*tstep);
			dTw3(i,dTwdt(i)*tstep);
			hi(i,horig(i)+dh3(i));
			Tt(i,Ttorig(i)+dTt3(i));
			Tw(i,Tworig(i)+dTw3(i));
		}
		fnStateDerivs(P,hi,Tt,Tw,mr0,Qr0,dPdt,dhdt,dTtdt,dTwdt,M1,VRHOU_1);
		dP4		=	dPdt*tstep;
		P		=	Porig + dP1/6.0 + dP2/3.0 + dP3/3.0 + dP4/6.0;
		M2	=	0.0;
		VRHOU_1	= 0.0;
		for(int i=1;i<=nnodes;i++)
		{
			dh4(i,dhdt(i)*tstep);
			dTt4(i,dTtdt(i)*tstep);
			dTw4(i,dTwdt(i)*tstep);
			hi(i,horig(i) + dh1(i)/6.0 + dh2(i)/3.0 + dh3(i)/3.0 + dh4(i)/6.0);
			Tt(i,Ttorig(i) + dTt1(i)/6.0 + dTt2(i)/3.0 + dTt3(i)/3.0 + dTt4(i)/6.0);
			Tw(i,Tworig(i) + dTw1(i)/6.0 + dTw2(i)/3.0 + dTw3(i)/3.0 + dTw4(i)/6.0);

			if(!state.setstate(P,hi(i)))
			{
				sprintf(errstr,"Cannot set properties at %fkPa and %fkJ/kg",P,hi(i));
				errorlog.Add("FVShellTubeHX::fnReEvalState ",errstr);
			}
			M2	+=	vol*state.getrho();
			VRHOU_1	+=	vol*state.getrho()*state.getu();
		}

		dDYNREFQTY	=	M2;
		dTOTMASSIN	+=	mrin*tstep;
		dTOTMASSOUT	+=	mrout*tstep;
		dREFMASSIMB	=	(dDYNREFQTY - dORIGREFQTY)-(dTOTMASSIN-dTOTMASSOUT);

		tstart	+=	tstep;
		Porig = P;
		horig = hi;
		Ttorig= Tt;
		Tworig= Tw;
	}//while tstart < tend

	//Energy Balance computation for the 1 sec step;
	double Pinit;
	Vector hinit;		hinit.make(nnodes);
	Vector Ttinit;		Ttinit.make(nnodes);
	Vector Twinit;		Twinit.make(nnodes);
	Pinit	=	VState(1);
	hinit	=	MState[1];
	Ttinit	=	MState[2];
	Twinit	=	MState[3];

	double tcap, wcap;
	tcap = VGeometry(9);
	wcap = VGeometry(10);

	double	Qrin, Qrout, Qwin, Qwout, Qrstored, Qtstored, Qwstored, Qtotstored;
	Cpw		=	4.1868;	//kJ/kg-K
	Qrin	=	mrin*hrin;		Qrout	=	mrout*hinit(nnodes);
	Qwin	=	mwat*Cpw*Twin;	Qwout	=	mwat*Cpw*Twinit(1);
	Qrstored=	VRHOU_1-VRHOU_0;
	Qtstored = 0.0;
	Qwstored = 0.0;
	for(int i=1;i<=nnodes;i++)
	{
		Qtstored += tcap*(Tt(i)-Ttinit(i));
		Qwstored += wcap*(Tw(i)-Twinit(i));
	}
	Qtotstored	=	Qrstored + Qtstored + Qwstored;
	HXEnBal	=	(Qrin + Qwin) - (Qrout+Qwout) - Qtotstored;
	//End of Energy Balance computation

	newstate1 = P;
	newstate2.putcol(1,hi);
	newstate2.putcol(2,Tt);
	newstate2.putcol(3,Tw);
	newstate2.putcol(4,mr0);
	newstate2.putcol(5,Qr0);
	
	return true;
}

double	FVShellTubeHX::fnRefMassBalance(double mgj_1, double rhogi, double rho, double tstep)
{
	//
	//Compute mass flow rate at jth interface
	//
	double vol;
	vol = VGeometry(6);

	return mgj_1 - (rhogi - rho)*vol/tstep;
}

double	FVShellTubeHX::fnRefEnergyBalance(double hgi_1, double hgi1, double Ttgi, double Trgi, double mgj, 
									   double mgj_1, double rho, double hi, double Hr, double tstep)
{
	//
	//Compute enthalpy of ith cell
	//
	double vol, htarea;
	vol = VGeometry(6);
	htarea = VGeometry(8);

	double upA, upB, upC;	//Upwind coefficients

	upA = (mgj_1>0)?mgj_1:0.0;
	upB = (-mgj>0)?(-mgj):0.0;
	upC = upA + upB + rho*vol/tstep;

	return (upA*hgi_1 + upB*hgi1 - Hr*htarea*(Trgi - Ttgi) + rho*vol*hi/tstep)/upC;
}

double	FVShellTubeHX::fnTubeEnergyBalance(double Hr, double Trg, double Hw, double Twg, double Tt,double tstep)
{
	double tubeid, tubeod;
	tubeid = VGeometry(1);
	tubeod = VGeometry(2);
	double htarea, ihtarea;
	htarea = VGeometry(8);
	ihtarea = htarea*tubeid/tubeod;

	double TubeCap;
	TubeCap	= VGeometry(9);

	double num, den;
	num = (Hr*htarea*Trg) + (Hw*ihtarea*Twg) + TubeCap*Tt/tstep;
	den = (Hr*htarea) + (Hw*ihtarea) + TubeCap/tstep;

	return num/den;
}

double	FVShellTubeHX::fnWaterEnergyBalance(double Twgi1, double Ttgi, double Tw, double Hw, double tstep)
{
	double tubeid, tubeod;
	tubeid = VGeometry(1);
	tubeod = VGeometry(2);
	double htarea, ihtarea;
	htarea = VGeometry(8);
	ihtarea = htarea*tubeid/tubeod;

	double mwat, Cpw;
	mwat = VINPUTS(1);
	Cpw  = 4.1868;	//kJ/kg-K
	double WaterCap;
	WaterCap = VGeometry(10);

	double num, den;
	num = (-mwat*Cpw)*Twgi1 + (Hw*ihtarea*Ttgi) + WaterCap*Tw/tstep;
	den = (-mwat*Cpw) + WaterCap/tstep + (Hw*ihtarea);

	return num/den;
}

bool	FVShellTubeHX::fnNodalMassEnergySolver(int nodeno, double Pg, double rho, Matrix& newstate2, double tstep)
{
	//
	//Set-up and evaluate mass and energy balance for cell specified by nodeno
	//
	int nnodes;
	nnodes = NNODES;

	char errmsg[512];

	//Capture relevant geometry
	double farea, htarea, tubeid,foulfactor;
	tubeid	=	VGeometry(1);
	farea	=	VGeometry(7);
	htarea	=	VGeometry(8);
	foulfactor	=	1 - VGeometry(4);

	//Capture relevant state info
	Vector hi;	hi.make(nnodes);	hi = MState[1];
	Vector Tt;	Tt.make(nnodes);	Tt = MState[2];
	Vector Tw;	Tw.make(nnodes);	Tw = MState[3];
	Vector Qr;	Qr.make(nnodes);	Qr = MState[5];

	//Capture boundary conditions
	double mwat, Twin, mrin, mrout,hrin;
	mwat = VINPUTS(1);	Twin = VINPUTS(2);	mrin = VINPUTS(3);	mrout = VINPUTS(4);	hrin = VINPUTS(5);

	r134astate state;
	double hgi, hgi_1, hgi1, Trgi, Hr, Hw,rhogg, mgj, mgj_1;
	double vel, Qflux, Qrg;
	double Ttgi,Twgi, Twgi1;
	int iters = 0;
	bool Converged = false;
	double URF = 1.0;

	//Current node guesses
	hgi		=	newstate2(nodeno,1);
	Ttgi	=	newstate2(nodeno,2);
	Twgi	=	newstate2(nodeno,3);
	mgj		=	newstate2(nodeno,4);
	Qrg		=	newstate2(nodeno,5);

	//Adjacent node states
	mgj_1	=	(nodeno==1)?mrin:newstate2(nodeno-1,4);
	hgi_1	=	(nodeno==1)?hrin:newstate2(nodeno-1,1);
	hgi1	=	(nodeno==nnodes)?newstate2(nnodes,1):newstate2(nodeno+1,1);
	Twgi1	=	(nodeno==nnodes)?Twin:newstate2(nodeno+1,3);

	if(!state.setstate(Pg,hgi))
	{
		sprintf(errmsg,"Could not set refrigerant state in node %d with %f kpa and %f kJ/kg",nodeno,Pg,hgi);
		errorlog.Add("FVShellTubeHX::fnNodalMassEnergySolver ",errmsg);
		return false;
	}
	rhogg	=	state.getrho();
	Trgi	=	state.getT();
	//Refrigerant mass-balance residual
	mgj		=	mgj - URF*(mgj - fnRefMassBalance(mgj_1,rhogg,rho,tstep));
	{
		//Compute refrigerant-side heat transfer coefficient
		vel		=	fabs(mgj+mgj_1)/(2*farea*rhogg);
		Qflux	=	fabs(Qrg/htarea);
		rpropertylist(1,state.getTs());		rpropertylist(3,state.getmu());		rpropertylist(4,state.getCp());
		rpropertylist(5,state.getk());		rpropertylist(6,state.getrho());	rpropertylist(7,state.getphase());
		rpropertylist(8,state.gethv());		rpropertylist(9,state.gethl());		rpropertylist(10,state.getCpl());
		rpropertylist(11,state.getrhol());	rpropertylist(12,state.getrhov());	rpropertylist(13,state.getkl());
		rpropertylist(14,state.getmul());	rpropertylist(15,vel);				rpropertylist(16,Qflux);
		rpropertylist(17,state.getx());		rpropertylist(18,state.getmuv());	rpropertylist(19,state.getCpv());
		rpropertylist(20,state.getkv());
		Hr		=	fnRefHTCoeff(rpropertylist,newstate2(nodeno,2));
	}
	//Refrigerant energy-balance residual
	hgi		=	hgi		-	URF*(hgi - fnRefEnergyBalance(hgi_1,hgi1,Ttgi,Trgi,mgj,mgj_1,rho,hi(nodeno),Hr,tstep));
	{
		//Compute water-side heat transfer coefficient
		watstate water;
		double area, vel;
		area	= MYPI*tubeid*tubeid/4.0;
		water.setstate(Twgi);
		vel		=	2.0*(mwat/NTUBES)/(water.getrho()*area);
		wpropertylist(1,water.getmu());
		wpropertylist(2,water.getk());
		wpropertylist(3,water.getCp());
		wpropertylist(4,water.getrho());
		wpropertylist(5,vel);
		Hw = foulfactor*fnWatHTCoeff(wpropertylist);
	}
	//Tube material energy-balance residual
	Ttgi	=	Ttgi	-	URF*(Ttgi - fnTubeEnergyBalance(Hr,Trgi,Hw,Twgi,Tt(nodeno),tstep));
	//Water energy-balance residual
	Twgi	=	Twgi	-	URF*(Twgi - fnWaterEnergyBalance(Twgi1,Ttgi,Tw(nodeno),Hw,tstep));

	Qrg		=	Hr*htarea*(Trgi - Ttgi);

	newstate2(nodeno,1,hgi);
	newstate2(nodeno,2,Ttgi);
	newstate2(nodeno,3,Twgi);
	newstate2(nodeno,4,mgj);
	newstate2(nodeno,5,Qrg);

	return true;
}

bool	FVShellTubeHX::fnHXMassEnergyBalancesSweep(double Pg, Vector& rho, Matrix& newstate2, double tstep)
{
	//
	//Perform one sweep of mass energy balances along heat-exchanger
	//
	int nnodes;
	nnodes = NNODES;

	char errmsg[512];

	for(int i=1;i<=nnodes;i++)
	{
		if(!fnNodalMassEnergySolver(i,Pg,rho(i),newstate2,tstep))
		{
			sprintf(errmsg,"Nodal mass-energy solver failed in node %d",i);
			errorlog.Add("FVShellTubeHX::fnHXMassEnergyBalancesSweep ",errmsg);
			return false;
		}
	}

	return true;
}

bool	FVShellTubeHX::fnEnthalpyProfileIterator(double Pg, Vector& rho, Matrix& newstate2, double tstep)
{
	//
	//Iterate on enthalpy and mass flow-rate profiles and check for convergence
	//
	int nnodes;
	nnodes = NNODES;

	//Capture state
	Vector hi;	hi.make(nnodes);	hi = MState[1];
	Vector Tt;	Tt.make(nnodes);	Tt = MState[2];
	Vector Tw;	Tw.make(nnodes);	Tw = MState[3];
	Vector mr;	mr.make(nnodes);	mr = MState[4];
	Vector Qr;	Qr.make(nnodes);	Qr = MState[5];

	bool Converged;
	Converged = false;
	Vector hg;	hg.make(nnodes);	hg		=	hi;//		hg.scale(0.95);
	Vector mrg;	mrg.make(nnodes);	mrg		=	mr;//		mrg.scale(1.05);
	Vector Ttg;	Ttg.make(nnodes);	Ttg		=	Tt;//		Ttg.scale(1.05);
	Vector Twg;	Twg.make(nnodes);	Twg		=	Tw;//		Twg.scale(1.05);

	newstate2.putcol(1,hg);
	newstate2.putcol(2,Ttg);
	newstate2.putcol(3,Twg);
	newstate2.putcol(4,mrg);
	newstate2.putcol(5,Qr);

	Vector hChange;	hChange.make(nnodes);
	double hgui;

	int iters=0;

	while(!Converged)
	{
		iters++;
		if(!fnHXMassEnergyBalancesSweep(Pg,rho,newstate2,tstep))
		{
			errorlog.Add("FVShellTubeHX::fnEnthalpyProfileIterator ","Mass-Energy balances sweep could not be completed succesfully.");
			return false;
		}
		for(int i=1;i<=nnodes;i++)
		{
			hgui = newstate2(i,1);
			hChange(i,fabs(hg(i)-hgui));
		}
		if(hChange.max()<=HTOL)
			Converged = true;
		else
		{
			//Save current profiles for comparison with next iteration
			newstate2.getcol(1,hg);
			newstate2.getcol(4,mrg);
		}
		if(iters>MAXITERS)
		{
			errorlog.Add("FVShellTubeHX::fnEnthalpyProfileIterator ","Iteration limit reached without convergence.");
			return false;
		}
	}//while !Converged

	return true;
}

bool	FVShellTubeHX::fnPressureIterator(double& newstate1, Matrix& newstate2, double& HXEnBal, double tstep)
{
	//
	//Iterate on pressure to achieve specified exit flow-rate
	//
	//Capture relevant geometry
	int nnodes;
	nnodes = NNODES;

	char errmsg[512];

	//Capture relevant boundary condition
	double mrout;
	mrout = VINPUTS(4);

	//Capture state
	double P;
	P	=	VState(1);
	Vector hi;		hi.make(nnodes);	hi	=	MState[1];
	Vector rho;		rho.make(nnodes);
	r134astate	state;
	for(int i=1;i<=nnodes;i++)
	{
		if(!state.setstate(P,hi(i)))
		{
			errorlog.Add("FVShellTubeHX::fnPressureIterator ","Could not set refrigerant state.");
			return false;
		}
		rho(i,state.getrho());
	}

	double Pg1, Pg2;
	Pg1 = P + 1.0;
	Pg2 = P + 2.0;

	double res1, res2, mroutg1, mroutg2;
	int iters=0;
	bool Converged = false;

	if(!fnEnthalpyProfileIterator(Pg1,rho,newstate2,tstep))
	{
		errorlog.Add("FVShellTubeHX::fnPressureIterator ","Enthalpy profile iterator failed in first secant pass.");
		return false;
	}
	mroutg1 = newstate2(nnodes,4);
	res1 = mrout - mroutg1;
	while(!Converged)
	{
		iters++;
		if(!fnEnthalpyProfileIterator(Pg2,rho,newstate2,tstep))
		{
			sprintf(errmsg,"Enthalpy profile iterator failed in 2nd secant pass at %d iteration.",iters);
			errorlog.Add("FVShellTubeHX::fnPressureIterator ",errmsg);
			return false;
		}
		mroutg2 = newstate2(nnodes,4);
		res2 = mrout - mroutg2;
		if(fabs(res1)<0.01&&fabs(res2)<0.01)
			Converged = true;
		else
		{
			double temp;
			//Updation of guess values
			temp	=	Pg1 - res1*(Pg1 - Pg2)/(res1 - res2);
			Pg1		=	Pg2;
			Pg2		=	temp;
			res1	=	res2;
		}
		if(iters>MAXITERS)
		{
			errorlog.Add("FVShellTubeHX::fnPressureIterator ","Iteration limit reached without convergence.");
			return false;
		}
	}//while !Converged
	newstate1 = (Pg1 + Pg2)/2.0;

	return true;
}

bool	FVShellTubeHX::fnReEvalStateImpE(double& newstate1, Matrix& newstate2, double& HXEnBal)
{

	double time=0, tstart=0.0, tend=1.0, tstep;
	char errmsg[512];

	while(tstart<tend)
	{
//		tstep	=	((tstart+dTSTEP)>tend)?(tend-tstart):dTSTEP;
		tstep = 0.5;
		if(!fnPressureIterator(newstate1,newstate2,HXEnBal,tstep))
		{
			sprintf(errmsg,"Pressure iteration failed at time %f",tstart);
			errorlog.Add("FVShellTubeHX::fnReEvalStateImpE ",errmsg);
			return false;
		}
		//Update state
		VState(1,newstate1);
		for(int i=1;i<=iSCOLS;i++)
			newstate2.getcol(i,MState[i]);
		tstart += tstep;
	}

	return true;
}

bool	FVShellTubeHX::fnAdvance1sec(int METHOD)
{
	char errstr[512];
	dTSTEP = (iHXTYPE==COND)?TCSTEP:TESTEP;
	dTSTEP	=	(dTIME<SUTIME)?DEFDT:dTSTEP;
	bTSIZING	=	false;

	Matrix	newstate2;
	newstate2.make(NNODES,iSCOLS);

	//Evaluate state after 1 second
	double Pnew;
	double HXEnBal;
	switch(METHOD)
	{
	case EXPE:
		if(!fnReEvalStateExpE(Pnew,newstate2,HXEnBal))		return false;
		break;
	case IMPEPC:
		if(!fnReEvalStateImpEPC(Pnew,newstate2,HXEnBal))	return false;
		break;
	case EXPEPC:
		if(!fnReEvalStateExpEPC(Pnew,newstate2,HXEnBal))	return false;
		break;
	case RK4:
		if(!fnReEvalStateRK4(Pnew,newstate2,HXEnBal))		return false;
		break;
	case IMPE:
		if(!fnReEvalStateImpE(Pnew,newstate2,HXEnBal))		return false;
		break;
	default:
		errorlog.Add("FVShellTubeHX::fnAdvance1sec ","Unknown method code.");
		break;
	}
	dTIME	+=	1.0;

	//Update state...part 1
	for(int i=1;i<=iSCOLS;i++)
		newstate2.getcol(i,MState[i]);

	//Leaving water temperature
	double Two;
	Two = MState[3](1);

	//Overall heat transfer
	double mwat,Twin,Cpw=4.1868,Qhx;
	mwat	=	VINPUTS(1);
	Twin	=	VINPUTS(2);
	Qhx		=	mwat*Cpw*(Two - Twin);

	//Exit enthalpy
	r134astate	exitstate;
	double P, hexit, Texit, Ts;
	int	phase;
	P		=	Pnew;
	hexit	=	MState[1](NNODES);
	if(!exitstate.setstate(P,hexit))
	{
		sprintf(errstr,"Could not set refrigerant properties with %f kPa and %f kJ/kg",P,hexit);
		errorlog.Add("FVShellTubeHX::fnAdvance1sec ",errstr);
		return false;
	}
	phase	=	exitstate.getphase();
	Ts		=	exitstate.getTs();
	Texit	=	exitstate.getT();
	double deltaT;
	switch(iHXTYPE)
	{
	case EVAP:
		deltaT	=	Texit-Ts;
		break;
	case COND:
		deltaT	=	Ts - Texit;
		break;
	}//switch

	//Update state...part 2
	VState(1,Pnew);		//Pressure
	VState(2,hexit);	//Enthalpy
	VState(3,Two);		//Leaving water temperature

    VOUTPUTS(1,Pnew);			VOUTPUTS(2,hexit);			VOUTPUTS(3,Two);
	VOUTPUTS(4,Qhx);			VOUTPUTS(5,deltaT);			VOUTPUTS(6,HXEnBal);
	VOUTPUTS(7,dDYNREFQTY);		VOUTPUTS(8,dREFMASSIMB);	VOUTPUTS(9,MState[1](BULBLOCN));

	return true;

}
//*******************************************************************************************

MicroTech::MicroTech()
{
	VGeometry.make(8);
	VState.make(1);
	iACTIONMODE = READY;
}

void	MicroTech::fnDefineGeometry(Vector& geom)
{
	
	double fopent, fcloset, startlim, rutime;
	fopent	= geom(3);
	fcloset	= geom(4);
	startlim= geom(5);
	rutime	= geom(6);

	VGeometry(1,geom(1));					//modlimit
	VGeometry(2,geom(2));					//deadband
	VGeometry(3,1/fopent);					//upmax
	VGeometry(4,1/fcloset);					//dnmax
	VGeometry(5,(100.0-startlim)/rutime);	//ruslope
	VGeometry(6,startlim);					//startlim
	VGeometry(7,geom(7));					//maxstep
	VGeometry(8,geom(8));					//waittime

}

void	MicroTech::fnSetState(double state, int mode)
{
	iRUNMODE	=	mode;	//STARTUP, NONSTARTUP
	dTIMER	=	0.0;
	VState(1,state);	//gamma
}

void	MicroTech::fnGetState(double& gamma)
{
	gamma	= VState(1);
}

double	MicroTech::fnReEvalState(Vector& bcs)
{
	double Tewoset,Tewo,RLApc,time;
	//Capture boundary conditions
	Tewoset = bcs(1);		//Set point chiller water temperature
	Tewo	= bcs(2);		//Actual chiller water temperature
	RLApc	= bcs(3);		//Rated Load Amps percent
	time	= bcs(4);		//time from startup

	double gamma;
	double modlimit, deadband,upmax,dnmax,maxstep,waittime;
	double ruslope,startlim;

	//Capture geometry
	modlimit = VGeometry(1);	deadband = VGeometry(2);	upmax	 = VGeometry(3);	dnmax	 = VGeometry(4);
	ruslope	 = VGeometry(5);	startlim = VGeometry(6);	maxstep	 = VGeometry(7);	waittime = VGeometry(8);

	//Capture state
	gamma	=	VState(1);
	double	GAMMA_MAX = 1.0, GAMMA_MIN = 0.05;

	double error, abserr;
	error	= Tewo - Tewoset;
	abserr	= fabs(error);

	if(iACTIONMODE==READY)
	{
		if(error>0)		{	bOPEN = true;	bCLOSE = false;	}
		else			{	bOPEN = false;	bCLOSE = true;	}
		if(abserr<=deadband)
			dACTIONTIME	= 0.0;
		else if(abserr>=deadband && abserr<=modlimit)
			dACTIONTIME	= ((abserr-deadband)/(modlimit-deadband))*maxstep;
		else
			dACTIONTIME	= waittime;

		iACTIONMODE = WAIT;
		dTIMER = 0.0;
	}
	if(iACTIONMODE==WAIT)
	{
		if(dACTIONTIME<1.0&&dACTIONTIME>0.0)
		{
			if((bOPEN)&&gamma<GAMMA_MAX)		gamma	+=	dACTIONTIME*upmax;
			if((bCLOSE)&&gamma>GAMMA_MIN)	gamma	-=	dACTIONTIME*dnmax;
			dACTIONTIME	=	0.0;
		}
		if(dACTIONTIME>=1.0)
		{
			if((bOPEN)&&gamma<GAMMA_MAX)		gamma	+=	upmax;
			if((bCLOSE)&&gamma>GAMMA_MIN)	gamma	-=	dnmax;
			dACTIONTIME	-=	1.0;
		}
		dTIMER	+=	1.0;
		if(dTIMER>=waittime)
		{
			iACTIONMODE = READY;
			dTIMER	=	0.0;
		}
	}
	return gamma;
}

bool	MicroTech::fnAdvance1sec(Vector& bcs, double& gamma)
{
	gamma = fnReEvalState(bcs);
	VState(1,gamma);
	return true;
}

//*******************************************************************************************

CentComp::CentComp()
{
	VGeometry.make(3);
	MGeometry = new Vector [4];
	VState.make(4);
	iN_INPUTS = 5;
	iN_OUTPUTS= 4;
	VINPUTS.make(iN_INPUTS);
	VOUTPUTS.make(iN_OUTPUTS);
}

bool	CentComp::fnDefineGeometry(const char* geomfile)
{
	ifstream	gfile(geomfile);
	if(!gfile)	return false;

	double etam, modlimit, deadband,fopent,fcloset,startlim,rutime,maxstep,waittime;
	gfile >> etam >> modlimit >> deadband >> fopent >> fcloset >> startlim >> rutime >> maxstep >> waittime;
	if(etam<=0||etam>1.0||modlimit<=0||modlimit<=deadband||deadband<=0||fopent<=0||fcloset<=0||
		startlim<=0||startlim>100||rutime<=0)
	{
		errorlog.Add("CentComp::fnDefineGeometry ","Invalid parameter values found in geometry file.");
		return false;
	}
	Vector contgeom;
	contgeom.make(8);
	contgeom(1,modlimit);	contgeom(2,deadband);	contgeom(3,fopent);	contgeom(4,fcloset);
	contgeom(5,startlim);	contgeom(6,rutime);		contgeom(7,maxstep);contgeom(8,waittime);

	Controller.fnDefineGeometry(contgeom);

	VGeometry(1,etam);		//Electro-mechanical efficiency
	VGeometry(2,rutime);	//Ramp up time
	VGeometry(3,startlim);	//Starting RLA limit
	MGeometry[1].make(5);	//max capacity map
	MGeometry[2].make(5);	//polytropic efficiency map
	MGeometry[3].make(3);	//RLA map

	double a1,a2,a3,a4,a5;
	//Max capacity map
	gfile >> a1 >> a2 >> a3 >> a4 >> a5;
	if(gfile.eof())	
	{	
		gfile.close();	
		errorlog.Add("CentComp::fnDefineGeometry ","End of file detected before construction could be completed.");
		return false;
	}
	MGeometry[1](1,a1);	MGeometry[1](2,a2);	MGeometry[1](3,a3);	MGeometry[1](4,a4);	MGeometry[1](5,a5);
	//Polytropic efficiency map
	gfile >> a1 >> a2 >> a3 >> a4 >> a5;
	if(gfile.eof())	
	{	
		gfile.close();	
		errorlog.Add("CentComp::fnDefineGeometry ","End of file detected before construction could be completed.");
		return false;
	}
	MGeometry[2](1,a1);	MGeometry[2](2,a2);	MGeometry[2](3,a3);	MGeometry[2](4,a4);	MGeometry[2](5,a5);
	//RLA map
	gfile >> a1 >> a2 >> a3;
	MGeometry[3](1,a1);	MGeometry[3](2,a2);	MGeometry[3](3,a3);

	gfile.close();
	return true;
}

bool	CentComp::fnSetState(double time, int mode, Vector& svs)
{
	dTIME	=	time;
	iMODE	=	mode;

	if(svs.size()!=5)
	{
		cerr << "Compressor initialization vector is size " << svs.size() << " instead of 5" << endl;
		return false;
	}
	VState(1,svs(1));	//Exit enthalpy
	VState(2,svs(2));	//Flow rate
	VState(3,svs(3));	//Power
	VState(4,svs(4));	//Losses
	Controller.fnSetState(svs(5),iMODE);

	return true;
}

bool	CentComp::fnSetState(double h2)
{
	dTIME	=	0.0;
	iMODE	=	STARTUP;

	VState(1,h2);
	VState(2,0.0);
	VState(3,0.0);
	VState(4,0.0);
	Controller.fnSetState(0.05,iMODE);

	return true;
}

void	CentComp::fnGetState(int& mode, Vector&	state)
{
	mode = iMODE;

	state.make(5);

	state(1,VState(1));	//Exit enthalpy
	state(2,VState(2));	//Flow rate
	state(3,VState(3));	//Power
	state(4,VState(4));	//Losses

	double gamma;
	Controller.fnGetState(gamma);
	state(5,gamma);			//Gamma
}

double	CentComp::fnWOVFlowRate(double P1, double T1, double P2)
{
	double	MCM[5];
	MCM[0] = MGeometry[1](1);
	MCM[1] = MGeometry[1](2);
	MCM[2] = MGeometry[1](3);
	MCM[3] = MGeometry[1](4);
	MCM[4] = MGeometry[1](5);

	return MCM[0] + MCM[1]*P1 + MCM[2]*P2 + MCM[3]*T1 + MCM[4]*P1*P2;
}

double	CentComp::fnPolytropicEfficiency(double Wpoly, double Vdotr)
{
	double PEFF[5], etap;
	PEFF[0] = MGeometry[2](1);
	PEFF[1] = MGeometry[2](2);
	PEFF[2] = MGeometry[2](3);
	PEFF[3] = MGeometry[2](4);	
	PEFF[4] = MGeometry[2](5);

	etap	= PEFF[0] + PEFF[1]*Vdotr + PEFF[2]*Vdotr*Vdotr + PEFF[3]*Wpoly + PEFF[4]*Wpoly*Wpoly;
	etap	= (etap>1.0)?1.0:etap;
	return (etap<ETAPMIN)?ETAPMIN:etap;
}

double	CentComp::fnRatedLoadAmpsPerCent(double Pwr)
{
	double RLA[3];
	RLA[0] = MGeometry[3](1);
	RLA[1] = MGeometry[3](2);
	RLA[2] = MGeometry[3](3);

	return	RLA[0] + RLA[1]*Pwr + RLA[2]*Pwr*Pwr;
}

double	CentComp::fnPolytropicWork(double P1, double v1, double P2, double v2)
{
	if(P2!=P1)
		return (P2*v2 - P1*v1)*log(P2/P1)/log((P2*v2)/(P1*v1));
	else
		return 0.0;
}
double	CentComp::fnInvRatedLoadAmpsPerCent(double RLA)
{
	double RLAM[3];
	RLAM[0] = MGeometry[3](1);
	RLAM[1] = MGeometry[3](2);
	RLAM[2] = MGeometry[3](3);

	double Pguess, disc;

	disc = sqrt(RLAM[1]*RLAM[1] - 4.0*RLAM[2]*(RLAM[0]-RLA));

	Pguess	=	(disc - RLAM[1])/(2.0*RLAM[2]);

	return Pguess;
}

bool	CentComp::fnResidue1(double P1, double h1, double P2, double h2, double mdotr, double& h2res)
{
	char errstr[512];
	r134astate outlet, inlet;
	if(!outlet.setstate(P2,h2))
	{
		sprintf(errstr,"Could not set refrigerant properties with %fkPa and %fkJ/kg",P2,h2);
		errorlog.Add("CentComp::fnRresidue1 ",errstr);
		return false;
	}
	if(!inlet.setstate(P1,h1))
	{
		sprintf(errstr,"Could not set refrigerant properties with %fkPa and %fkJ/kg",P1,h1);
		errorlog.Add("CentComp::fnRresidue1 ",errstr);
		return false;
	}

	double v1, v2, Wp, Vdotr, etap, Wc, h2new;
	v1		=	inlet.getv();
	v2		=	outlet.getv();
	Wp		=	fnPolytropicWork(P1,v1,P2,v2);
	Vdotr	=	mdotr*v1;
	etap	=	fnPolytropicEfficiency(Wp,Vdotr);
	Wc		=	Wp/etap;
	h2new	=	Wc + h1;
	h2res	=	h2 - h2new;

	return true;
}

bool	CentComp::fnResidue2(double Pwr, double h2, double mdotr, double h1, double& res)
{
	double etam;
	etam = VGeometry(1);
	res	= Pwr - mdotr*(h2-h1)/etam;
	return true;
}
/*
bool	CentComp::fnQuasiSteadyState1D(double P1, double h1, double P2, double mdotr, double h2guess, double& h2)
{
	char errstr[512];
	double h2g1, h2g2, h2res1, h2res2, temp;
	int iters=0;
	bool	Converged = false;

	h2g1	=	h2guess-2.0;
	h2g2	=	h2guess;

	if(P2!=P1)
	{
		//With guess value 1
		if(!fnResidue1(P1,h1,P2,h2g1,mdotr,h2res1))
		{
			sprintf(errstr,"Residue1 could not be computed with P1=%fkPa,h1=%fkJ/kg,P2=%fkPa,h2g1=%fkJ/kg,mc=%fkg/s.\nResidual remains %f",P1,h1,P2,h2g1,mdotr,h2res1);
			errorlog.Add("CentComp::fnQuasiSteadyState1D ",errstr);
			return false;
		}
		if(fabs(h2res1)<=HTOL)
		{
			h2 = h2g1;
			return true;
		}
		while(!Converged)
		{	//Iterate for h2 using Secant method
			//With guess value 2
			if(!fnResidue1(P1,h1,P2,h2g2,mdotr,h2res2))
			{
				sprintf(errstr,"Residue1 could not be computed with P1=%fkPa,h1=%fkJ/kg,P2=%fkPa,h2g2=%fkJ/kg,mc=%fkg/s.\nResidual remains %f",P1,h1,P2,h2g2,mdotr,h2res2);
				errorlog.Add("CentComp::fnQuasiSteadyState1D ",errstr);
				return false;
			}
			if(fabs(h2res2)<=HTOL)
				Converged = true;
			else
			{
				//Updation of guess values
				temp	=	h2g1 - h2res1*(h2g1 - h2g2)/(h2res1 - h2res2);
				h2g1	=	h2g2;
				h2g2	=	temp;
				h2res1	=	h2res2;
				if(iters++>MAXITERS)
				{
					sprintf(errstr,"Convergence not achieved in %d iterations.",MAXITERS);
					errorlog.Add("CentComp::fnQuasiSteadyState1D ",errstr);
					return false;
				}
			}
		} //test for convergence
		h2	=	(h2g1+h2g2)/2;
	}
	else
		h2 = h1;

	return true;
}
*/
bool	CentComp::fnQuasiSteadyState1D(double P1, double h1, double P2, double mdotr, double h2guess, double& h2)
{
	char errstr[512];
	double h2g1, h2g2, h2gmid, h2res1, h2res2, h2resmid;
	int iters=0;
	bool Bisected = false;
	bool Bracketed = false;

	h2g1	=	h2guess;
	h2g2	=	h2guess+10.0;

	if(P2!=P1)
	{
	//Solve for h2 using Bracketing and Bisection
		//With guess value 1
		while(!Bracketed)
		{
			if(!fnResidue1(P1,h1,P2,h2g1,mdotr,h2res1))
			{
				sprintf(errstr,"Bracketing failed at 1st limit. P1=%f, h1=%f, P2=%f, h2g1=%f, mdotr=%f, h2res1=%f",P1,h1,P2,h2g1,mdotr,h2res1);
				errorlog.Add("CentComp::fnQuasiSteadyState1D ",errstr);
				return false;
			}
			if(!fnResidue1(P1,h1,P2,h2g2,mdotr,h2res2))
			{
				sprintf(errstr,"Bracketing failed at 2nd limit. P1=%f, h1=%f, P2=%f, h2g2=%f, mdotr=%f, h2res2=%f",P1,h1,P2,h2g2,mdotr,h2res2);
				errorlog.Add("CentComp::fnQuasiSteadyState1D ",errstr);
				return false;
			}
			if(!SAMESIGN(h2res1,h2res2))
				Bracketed = true;
			else
			{
				if(fabs(h2res1) < fabs(h2res2))
					h2g1	+=	1.6*(h2g2-h2g1);
				else
					h2g2	+=	1.6*(h2g2-h2g1);
			}
		}//while !Bracketed
		double dx, rtb;
		rtb = (h2res1<0.0)?(dx=h2g2-h2g1,h2g1):(dx=h2g1-h2g2,h2g2);
		while(!Bisected)
		{
			dx	*=	0.5;
			h2gmid = rtb + dx;
			if(!fnResidue1(P1,h1,P2,h2gmid,mdotr,h2resmid))
			{
				sprintf(errstr,"Bracketing failed at 2nd limit. P1=%f, h1=%f, P2=%f, h2gmid=%f, mdotr=%f, h2resmid=%f",P1,h1,P2,h2gmid,mdotr,h2resmid);
				errorlog.Add("CentComp::fnQuasiSteadyState1D ",errstr);
				return false;
			}
			if(h2resmid<=0.0)
				rtb	=	h2gmid;
			if(fabs(dx)<HTOL || h2resmid ==0.0)
			{
				h2 = rtb;
				Bisected = true;
			}
		}
	}
	else
		h2 = h1;

	return true;
}

bool	CentComp::fnQuasiSteadyState2D(double P1, double h1, double P2, double Pwr, 
									   double h2g, double mdotrg, double& h2, double& mdotr)
{
	char errstr[512];
	Vector F;	F.make(2);		//Residual function vector
	Matrix J;	J.make(2,2);	//Jacobian
	Vector X;	X.make(2);		//Argument values
	Vector DX;	DX.make(2);		//Corrections
	Vector X1;	X1.make(2);		//Tentative new arguments
	Vector F1;	F1.make(2);		//Tentative new residual

	double res;
	int iters=0;

	X(1,h2g);	X(2,mdotrg);
	do
	{
		if(!fnResidue1(P1,h1,P2,X(1),X(2),res))	
		{
			sprintf(errstr,"Could not compute residue1 with P1=%fkPa,h1=%fkJ/kg,P2=%fkPa,h2g=%fkJ/kg,mc=%fkg/s.\nResidue remains %f",P1,h1,P2,X(1),X(2),res);
			errorlog.Add("CentComp::fnQuasiSteadyState2D ",errstr);
			return false;
		}
		F(1,res);
		if(!fnResidue2(Pwr,X(1),X(2),h1,res))
			return false;
		F(2,res);
		{
			if(!fnResidue1(P1,h1,P2,X(1)+DELTA,X(2),res))
			{
				sprintf(errstr,"Could not compute residue1 with P1=%fkPa,h1=%fkJ/kg,P2=%fkPa,h2g=%fkJ/kg,mc=%fkg/s.\nResidue remains %f",P1,h1,P2,X(1),X(2),res);
				errorlog.Add("CentComp::fnQuasiSteadyState2D ",errstr);
				return false;
			}
			res	= (res-F(1))/DELTA;
			J(1,1,res);
			if(!fnResidue1(P1,h1,P2,X(1),X(2)+DELTA,res))
			{
				sprintf(errstr,"Could not compute residue1 with P1=%fkPa,h1=%fkJ/kg,P2=%fkPa,h2g=%fkJ/kg,mc=%fkg/s.\nResidue remains %f",P1,h1,P2,X(1),X(2),res);
				errorlog.Add("CentComp::fnQuasiSteadyState2D ",errstr);
				return false;
			}
			res	= (res-F(1))/DELTA;
			J(1,2,res);
			if(!fnResidue2(Pwr,X(1)+DELTA,X(2),h1,res))		return false;
			res	= (res-F(2))/DELTA;
			J(2,1,res);
			if(!fnResidue2(Pwr,X(1),X(2)+DELTA,h1,res))		return false;
			res	= (res-F(2))/DELTA;
			J(2,2,res);
		}
		solveLU(J,F,DX);
		bool done = false;
		do{
			X1(1,X(1)-DX(1));
			X1(2,X(2)-DX(2));
			if(!fnResidue1(P1,h1,P2,X1(1),X1(2),res))
			{
				sprintf(errstr,"Could not compute residue1 with P1=%fkPa,h1=%fkJ/kg,P2=%fkPa,h2g=%fkJ/kg,mc=%fkg/s.\nResidue remains %f",P1,h1,P2,X(1),X(2),res);
				errorlog.Add("CentComp::fnQuasiSteadyState2D ",errstr);
				return false;
			}
			F1(1,res);
			if(!fnResidue2(Pwr,X1(1),X1(2),h1,res))		return false;
			F1(2,res);
			if(F1.norm()>F.norm())
			{
				DX(1,DX(1)*0.5);
				DX(2,DX(2)*0.5);
			}
			else 
			{
				X = X1;	F = F1;
				done = true;
			}
		}while(!done);
		iters++;
		if(iters>=MAXITERS)
		{
			sprintf(errstr,"Convergence not achieved in %d iterations",MAXITERS);
			errorlog.Add("CentComp::fnQuasiSteadyState2D ",errstr);
			return false;
		}
	}while(F.norm()>FTOL);

	h2		=	X(1);
	mdotr	=	X(2);

	return true;
}

bool	CentComp::fnReEvalState(Vector& newstate)
{
	char errstr[512];
	//Capture boundary conditions
	double P1, h1, P2, Tewo, TewoSet, time;
	P1		= VINPUTS(1);	h1		= VINPUTS(2);	P2		= VINPUTS(3);
	TewoSet	= VINPUTS(4);	Tewo	= VINPUTS(5);	time	= dTIME;

	//Run controller
	Vector	contbcs;	contbcs.make(4);
	contbcs(1,TewoSet);	contbcs(2,Tewo);
	double gamma;
	Controller.fnAdvance1sec(contbcs,gamma);

	//Calculate flow rate
	double mdotr, mrmax, T1, v1;
	r134astate inlet;
	if(!inlet.setstate(P1,h1))
	{
		sprintf(errstr,"Could not set refrigerant properties with %fkPa and %fkJ/kg",P1,h1);
		errorlog.Add("CentComp::fnReEvalState ",errstr);
		return false;
	}
	v1	=	inlet.getv();
	T1	=	inlet.getT();
	mrmax	=	fnWOVFlowRate(P1,T1,P2);
	mdotr	=	gamma*mrmax;

	//Capture state
	double h2, h2g;
	h2		=	VState(1);
	h2g		=	h2;

	if(!fnQuasiSteadyState1D(P1,h1,P2,mdotr,h2g,h2))	return false;
	double	etam, Qloss, Pwr;
	etam	=	VGeometry(1);
	Pwr		=	mdotr*(h2-h1)/etam;
	Qloss	=	(1-etam)*Pwr;
	
	double RLApc, RLAlimit, rutime, startlim, ruslope;
	rutime	=	VGeometry(2);
	startlim=	VGeometry(3);
	ruslope	=	(100.0-startlim)/rutime;
	RLAlimit=	(iMODE==STARTUP)?(startlim+time*ruslope):100.0;
	RLAlimit=	(RLAlimit>100.0)?100.0:RLAlimit;
	RLAlimit=	RLAlimit + 10.0;

	RLApc	=	fnRatedLoadAmpsPerCent(Pwr);
	if(RLApc>RLAlimit)
	{
		double mdotrg;
		mdotrg	=	VState(2);
		Pwr	= fnInvRatedLoadAmpsPerCent(RLAlimit);
		RLApc	= RLAlimit;
		if(!fnQuasiSteadyState2D(P1,h1,P2,Pwr,h2g,mdotrg,h2,mdotr))	return false;
		Qloss	=	(1-etam)*Pwr;
	}

	newstate(1,h2);
	newstate(2,mdotr);
	newstate(3,Pwr);
	newstate(4,Qloss);

	return true;
}

bool	CentComp::fnAdvance1sec()
{
	Vector	newstate;
	newstate.make(4);
	if(iMODE!=SHUTDOWN)
	{
		double rutime;
		rutime = VGeometry(2);
		if(dTIME>rutime&&iMODE==STARTUP)
			iMODE = NORMAL;
		if(!fnReEvalState(newstate))	return false;
	}
	else
	{
		newstate(1,200.0);
		newstate(2,0.0);
		newstate(3,0.0);
		newstate(4,0.0);
	}
	dTIME	+=	1.0;

	VState(1,newstate(1));
	VState(2,newstate(2));
	VState(3,newstate(3));
	VState(4,newstate(4));

	VOUTPUTS(1,newstate(1));	//Exit enthalpy
	VOUTPUTS(2,newstate(2));	//Flow rate
	VOUTPUTS(3,newstate(3));	//Power
	VOUTPUTS(4,newstate(4));	//Loss

	return true;
}

void	CentComp::shutoff()
{
	iMODE = SHUTDOWN;
}

//*******************************************************************************************

LumpedCapacitance::LumpedCapacitance()
{
	VGeometry.make(1);	//Lumped capacitance
	VState.make(1);	//Tbulb
}

bool	LumpedCapacitance::fnDefineGeometry(double C)
{
	VGeometry(1,C);
	return true;
}

void	LumpedCapacitance::fnSetState(double Tbulb)
{
	VState(1,Tbulb);
}

void	LumpedCapacitance::fnGetState(double& Tbulb)
{
	Tbulb = VState(1);
}

double	LumpedCapacitance::fnReEvalState(double T1)
{
	double	dTbdt;
	double	ti=0.0, tf=1.0, dt = 0.02;

	//Capture geometry
	double C;
	C = VGeometry(1);

	//Capture state
	double Tb;
	Tb = VState(1);

	while(ti<=tf)
	{
		dTbdt = (T1-Tb)/C;
		Tb +=	dTbdt*dt;
		ti += dt;
	}

	return Tb;
}

bool	LumpedCapacitance::fnAdvance1sec(double T1, double& Tb)
{
	Tb = fnReEvalState(T1);
	VState(1,Tb);
	return true;
}

//*******************************************************************************************

ThermoStaticExValve::ThermoStaticExValve()
{
	VGeometry.make(8);
	VState.make(3);
	iN_INPUTS = 4;
	iN_OUTPUTS= 3;
	VINPUTS.make(iN_INPUTS);
	VOUTPUTS.make(iN_OUTPUTS);
}

bool	ThermoStaticExValve::fnDefineGeometry(const char* geomfile)
{
	ifstream	gfile(geomfile);
	if(!gfile)	return false;

	double Amax, theta, C_d, kspring, Cbulb, dPclose;
	gfile >> Amax >> theta >> C_d >> kspring >> Cbulb >> dPclose;
	gfile.close();

	double	maxlift, am0, am1, ro;
	ro		=	sqrt(Amax/MYPI);									//units are 'm'
	am0		=	2*MYPI*ro*tan(theta*MYPI/180.0);					//units are 'm'
	am1		=	-MYPI*tan(theta*MYPI/180.0)*tan(theta*MYPI/180.0);	//units are 'ND'
	maxlift =	ro/tan(theta*MYPI/180.0);							//units are 'm'

	VGeometry(1,Amax);
	VGeometry(2,C_d);
	VGeometry(3,kspring);
	VGeometry(4,Cbulb);
	VGeometry(5,am0);
	VGeometry(6,am1);
	VGeometry(7,maxlift);
	VGeometry(8,dPclose);

	return true;
}

bool	ThermoStaticExValve::fnSetState(double time, double mdotv, double Tbulb)
{
	char errstr[512];
	dTIME = time;

	double Pbulb;
	r500state	bstate;
	if(!bstate.setstateT(Tbulb))
	{
		sprintf(errstr,"Could not set R500 properties with %fC",Tbulb);
		errorlog.Add("ThermoStaticExValve::fnSetState ",errstr);
		return false;
	}
	Pbulb = bstate.getP();

	VState(1,mdotv);
	VState(2,Pbulb);
	VState(3,Tbulb);

	return true;
}

bool	ThermoStaticExValve::fnSetState(double Tbulb)
{
	char errstr[512];
	dTIME	=	0.0;

	double	Pbulb;
	r500state	bstate;
	if(!bstate.setstateT(Tbulb))
	{
		sprintf(errstr,"Could not set R500 properties with %fC",Tbulb);
		errorlog.Add("ThermoStaticExValve::fnSetState ",errstr);
		return false;
	}
	Pbulb = bstate.getP();

	VState(1,0.0);
	VState(2,Pbulb);
	VState(3,Tbulb);

	return true;
}

void	ThermoStaticExValve::fnGetState(Vector& state)
{
	state.make(3);

	state(1,VState(1));	//mdotv
	state(2,VState(3));	//Tbulb
	state(3,VState(4));	//lift
}

bool	ThermoStaticExValve::fnReEvalState1()
{
	char errstr[512];

	//Capture boundary conditions
	double P1, h1, P2, h3;
	P1	=	VINPUTS(1);
	h1	=	VINPUTS(2);
	P2	=	VINPUTS(3);
	h3	=	VINPUTS(4);

	//Capture geometry
	double Amax, Cd, kspring, am0, am1, maxlift,Cbulb, dPclose;
	Amax	= VGeometry(1);
	Cd		= VGeometry(2);
	kspring	= VGeometry(3);
	Cbulb	= VGeometry(4);
	am0		= VGeometry(5);
	am1		= VGeometry(6);
	maxlift	= VGeometry(7);
	dPclose	= VGeometry(8);


	//Capture state
	double mdotv, Pbulb,Tbulb;
	mdotv	= VState(1);
	Pbulb	= VState(2);
	Tbulb	= VState(3);

	//Re-calculate state variables
	r134astate	evapout;
	if(!evapout.setstate(P1,h1))
	{
		sprintf(errstr,"Could not set refrigerant properties with %fkPa and %fkJ/kg",P1,h1);
		errorlog.Add("ThermoStaticExValve::fnReEvalState ",errstr);
		return false;
	}
	double	T1;
	T1	=	evapout.getT();

	r500state	bulb;

	//Valve dynamics here...
	double ti=0.0, tf=1.0, dt=0.1, dTbdt;
	while(ti<=tf)
	{
		dTbdt	=	(T1-Tbulb)/Cbulb;
		Tbulb	+=	dTbdt*dt;
		ti		+=	dt;
	}

	double lift, Area, dPsh;
	if(!bulb.setstateT(Tbulb))
	{
		sprintf(errstr,"Could not set R500 properties with %fC",Tbulb);
		errorlog.Add("ThermoStaticExValve::fnSetState ",errstr);
		return false;
	}
	Pbulb	=	bulb.getP();
	dPsh	=	Pbulb - P1 - dPclose;
	dPsh	=	(dPsh<0)?0.0:dPsh;

	lift	=	(kspring*dPsh);					//units are 'm'
	lift	=	(lift>maxlift)?maxlift:lift;	//units are 'm'
	Area	=	(am0*lift + am1*lift*lift);		//units are 'm^2'
	Area	=	(Area>Amax)?Amax:Area;

	double v3;
	r134astate	condout;
	if(!condout.setstate(P2,h3))
	{
		sprintf(errstr,"Could not set refrigerant properties with %fkPa and %fkJ/kg",P2,h3);
		errorlog.Add("ThermoStaticExValve::fnReEvalState ",errstr);
		return false;
	}
	v3	=	condout.getv();

	double dP	=	(P2-P1)*1e3;
	if(dP>=0)
		mdotv	=	Cd*Area*sqrt(dP/v3);
	else
		mdotv	=	0.0;

	VOUTPUTS(1,mdotv);
	VOUTPUTS(2,Pbulb);
	VOUTPUTS(3,Tbulb);
	VOUTPUTS(4,lift);
	VOUTPUTS(5,Area);

	return true;
}

bool	ThermoStaticExValve::fnReEvalState2()
{
	char errstr[512];

	//Capture boundary conditions
	double P1, h1, P2, h3;
	P1	=	VINPUTS(1);
	h1	=	VINPUTS(2);
	P2	=	VINPUTS(3);
	h3	=	VINPUTS(4);

	//Capture geometry
	double Cbulb;
	Cbulb	= VGeometry(4);

	//Capture state
	double mdotv, Tbulb;
	Tbulb	= VState(3);

	//Re-calculate state variables
	r134astate	evapout;
	if(!evapout.setstate(P1,h1))
	{
		sprintf(errstr,"Could not set refrigerant properties with %fkPa and %fkJ/kg",P1,h1);
		errorlog.Add("ThermoStaticExValve::fnReEvalState ",errstr);
		return false;
	}
	double	T1, Tsh, Tcsat, Tesat, T3, Tsub;
	T1	=	evapout.getT();
	Tesat=	evapout.getTs();

	r500state	bulb;
	//Valve dynamics here...
	double ti=0.0, tf=1.0, dt=0.1, dTbdt;
	while(ti<=tf)
	{
		dTbdt	=	(T1-Tbulb)/Cbulb;
		Tbulb	+=	dTbdt*dt;
		ti		+=	dt;
	}

	Tsh = Tbulb - Tesat;

	r134astate	condout;
	if(!condout.setstate(P2,h3))
	{
		sprintf(errstr,"Could not set refrigerant properties with %fkPa and %fkJ/kg",P2,h3);
		errorlog.Add("ThermoStaticExValve::fnReEvalState ",errstr);
		return false;
	}
	T3	=	condout.getT();
	Tcsat=	condout.getTs();
	Tsub=	Tcsat - T3;

	double dP	=	(P2-P1);

	if(dP>=0)
	{
		mdotv	=	5.47961709E-01+8.66305971E-02*Tesat-9.36078079E-02*Tesat*Tesat-5.64327093E-01*Tsub+5.84252560E-03*dP+1.43693444E-01*Tesat*Tsub;
		mdotv	=	(mdotv>0.0)?mdotv:0.0;
	}
	else
		mdotv	=	0.0;

	VOUTPUTS(1,mdotv);
//	VOUTPUTS(2,Pbulb);
	VOUTPUTS(3,Tbulb);
//	VOUTPUTS(4,lift);
//	VOUTPUTS(5,Area);

	return true;
}

bool	ThermoStaticExValve::fnReEvalState3()
{
	char errstr[512];

	//Capture boundary conditions
	double P1, h1, P2, h3;
	P1	=	VINPUTS(1);
	h1	=	VINPUTS(2);
	P2	=	VINPUTS(3);
	h3	=	VINPUTS(4);

	//Capture geometry
	double Amax, Cd, kspring, am0, am1, maxlift,Cbulb, dPclose;
	Amax	= VGeometry(1);
	Cd		= VGeometry(2);
	kspring	= VGeometry(3);
	Cbulb	= VGeometry(4);
	am0		= VGeometry(5);
	am1		= VGeometry(6);
	maxlift	= VGeometry(7);
	dPclose	= VGeometry(8);

	//Capture state
	double mdotv, Pbulb,Tbulb;
	mdotv	= VState(1);
	Pbulb	= VState(2);
	Tbulb	= VState(3);

	//Re-calculate state variables
	r134astate	evapout;
	if(!evapout.setstate(P1,h1))
	{
		sprintf(errstr,"Could not set refrigerant properties with %fkPa and %fkJ/kg",P1,h1);
		errorlog.Add("ThermoStaticExValve::fnReEvalState ",errstr);
		return false;
	}
	double	T1, Tesat, Tsh, Tbdrv;
	T1		=	evapout.getT();
	Tesat	=	evapout.getTs();
	Tsh		=	T1 - Tesat;
	Tbdrv	=	0.7865*Tsh + Tesat;

	r500state	bulb;
	//Valve dynamics here...
	double ti=0.0, tf=1.0, dt=0.1, dTbdt;
	while(ti<=tf)
	{
		dTbdt	=	(Tbdrv-Tbulb)/Cbulb;
		Tbulb	+=	dTbdt*dt;
		ti		+=	dt;
	}

	double lift, Area, dPsh;
	if(!bulb.setstateT(Tbulb))
	{
		sprintf(errstr,"Could not set R500 properties with %fC",Tbulb);
		errorlog.Add("ThermoStaticExValve::fnSetState ",errstr);
		return false;
	}
	Pbulb	=	bulb.getP();
	dPsh	=	Pbulb - P1 - dPclose;
	dPsh	=	(dPsh<0)?0.0:dPsh;

	lift	=	(kspring*dPsh);					//units are 'm'
	lift	=	(lift>maxlift)?maxlift:lift;	//units are 'm'
	Area	=	(am0*lift + am1*lift*lift);		//units are 'm^2'
	Area	=	(Area>Amax)?Amax:Area;

	double v3;
	r134astate	condout;
	if(!condout.setstate(P2,h3))
	{
		sprintf(errstr,"Could not set refrigerant properties with %fkPa and %fkJ/kg",P2,h3);
		errorlog.Add("ThermoStaticExValve::fnReEvalState ",errstr);
		return false;
	}
	v3	=	condout.getv();

	double dP	=	(P2-P1)*1e3;
	if(dP>=0)
		mdotv	=	Cd*Area*sqrt(dP/v3);
	else
		mdotv	=	0.0;

	VOUTPUTS(1,mdotv);
	VOUTPUTS(2,Pbulb);
	VOUTPUTS(3,Tbulb);
	VOUTPUTS(4,lift);
	VOUTPUTS(5,Area);

	return true;
}

bool	ThermoStaticExValve::fnReEvalState4()
{
	char errstr[512];

	//Capture boundary conditions
	double P1, h1, P2, h3;
	P1	=	VINPUTS(1);
	h1	=	VINPUTS(2);
	P2	=	VINPUTS(3);
	h3	=	VINPUTS(4);

	//Capture geometry
	double Amax, Cd, kspring, am0, am1, maxlift,Cbulb, dPclose;
	Amax	= VGeometry(1);
	Cd		= VGeometry(2);
	kspring	= VGeometry(3);
	Cbulb	= VGeometry(4);
	am0		= VGeometry(5);
	am1		= VGeometry(6);
	maxlift	= VGeometry(7);
	dPclose	= VGeometry(8);

	//Capture state
	double mdotv, Pbulb,Tbulb;
	mdotv	= VState(1);
	Pbulb	= VState(2);
	Tbulb	= VState(3);

	//Re-calculate state variables
	r134astate	evapout;
	if(!evapout.setstate(P1,h1))
	{
		sprintf(errstr,"Could not set refrigerant properties with %fkPa and %fkJ/kg",P1,h1);
		errorlog.Add("ThermoStaticExValve::fnReEvalState ",errstr);
		return false;
	}
	double	T1, Tesat, Tsh, Tbdrv;
	T1		=	evapout.getT();
	Tesat	=	evapout.getTs();
	Tsh		=	T1 - Tesat;
	Tbdrv	=	0.7865*Tsh + Tesat;

	//Valve dynamics here...
	double ti=0.0, tf=1.0, dt=0.1, dTbdt;
	while(ti<=tf)
	{
		dTbdt	=	(Tbdrv-Tbulb)/Cbulb;
		Tbulb	+=	dTbdt*dt;
		ti		+=	dt;
	}

	double lift, Area, dP;
	lift	=	5e-5*Tbulb*Tbulb - 4e-5*Tbulb + 0.0453;
	lift	=	(lift>maxlift)?maxlift:lift;	//units are 'm'
	Area	=	(am0*lift + am1*lift*lift);		//units are 'm^2'
	Area	=	(Area>Amax)?Amax:Area;

	double v3;
	r134astate	condout;
	if(!condout.setstate(P2,h3))
	{
		sprintf(errstr,"Could not set refrigerant properties with %fkPa and %fkJ/kg",P2,h3);
		errorlog.Add("ThermoStaticExValve::fnReEvalState ",errstr);
		return false;
	}
	v3	=	condout.getv();

	dP	=	(P2 - P1);
	dP	*=	1e3;
	if(dP>=0)
		mdotv	=	Cd*Area*sqrt(dP/v3);
	else
		mdotv	=	0.0;

	VOUTPUTS(1,mdotv);
	VOUTPUTS(2,Pbulb);
	VOUTPUTS(3,Tbulb);
	VOUTPUTS(4,lift);
	VOUTPUTS(5,Area);

	return true;
}

bool	ThermoStaticExValve::fnReEvalState5(double Tewo)
{
	char errstr[512];

	//Capture boundary conditions
	double P1, h1, P2, h3;
	P1	=	VINPUTS(1);
	h1	=	VINPUTS(2);
	P2	=	VINPUTS(3);
	h3	=	VINPUTS(4);

	//Capture geometry
	double Amax, Cd, Cbulb;
	Amax	= VGeometry(1);
	Cd		= VGeometry(2);
	Cbulb	= VGeometry(4);

	//Capture state
	double mdotv, Pbulb,Tbulb;
	mdotv	= VState(1);
	Pbulb	= VState(2);
	Tbulb	= VState(3);

	//Re-calculate state variables
	r134astate	evapout;
	if(!evapout.setstate(P1,h1))
	{
		sprintf(errstr,"Could not set refrigerant properties with %fkPa and %fkJ/kg",P1,h1);
		errorlog.Add("ThermoStaticExValve::fnReEvalState ",errstr);
		return false;
	}
	double	T1, Tesat, Tsh, Tbdrv;
	T1		=	evapout.getT();
	Tesat	=	evapout.getTs();
	Tsh		=	T1 - Tesat;
	Tbdrv	=	0.7865*Tsh + Tesat;

	//Valve dynamics here...
	double ti=0.0, tf=1.0, dt=0.1, dTbdt;
	while(ti<=tf)
	{
		dTbdt	=	(Tbdrv-Tbulb)/Cbulb;
		Tbulb	+=	dTbdt*dt;
		ti		+=	dt;
	}

	double Area;
	Area	=	-1.54572380E-04
				+8.30280214E-07*P1
				+4.23745424E-08*P2
				+1.74806056E-05*Tbulb
				-3.80888745E-06*Tbulb*Tbulb
				-3.57530725E-05*Tewo
				+3.41217208E-06*Tewo*Tewo;

	Area	=	(Area<0.0)?0.0:Area;

	double v3, dP;
	r134astate	condout;
	if(!condout.setstate(P2,h3))
	{
		sprintf(errstr,"Could not set refrigerant properties with %fkPa and %fkJ/kg",P2,h3);
		errorlog.Add("ThermoStaticExValve::fnReEvalState ",errstr);
		return false;
	}
	v3	=	condout.getv();

	dP	=	(P2 - P1);
	dP	*=	1e3;
	if(dP>=0)
		mdotv	=	Cd*Area*sqrt(dP/v3);
	else
		mdotv	=	0.0;

	VOUTPUTS(1,mdotv);
	VOUTPUTS(2,Pbulb);
	VOUTPUTS(3,Tbulb);
//	VOUTPUTS(4,lift);
	VOUTPUTS(5,Area);

	return true;
}

bool	ThermoStaticExValve::fnAdvance1sec(double Tewo)
{

		if(!fnReEvalState1())
			return false;

	dTIME	+=	1.0;

	VState(1,VOUTPUTS(1));	//mdotv
	VState(2,VOUTPUTS(2));	//Pbulb
	VState(3,VOUTPUTS(3));	//Tbulb


	return true;
}

//*******************************************************************************************

Orifice::Orifice()
{
	VGeometry.make(2);
	VState.make(1);
	iN_INPUTS = 3;
	iN_OUTPUTS= 1;
	VINPUTS.make(iN_INPUTS);
	VOUTPUTS.make(iN_OUTPUTS);
}

bool	Orifice::fnDefineGeometry(const char* geomfile)
{
	ifstream	gfile(geomfile);
	if(!gfile)	return false;

	double A, C_d;
	gfile >> A >> C_d;
	gfile.close();

	VGeometry(1,A);
	VGeometry(2,C_d);

	return true;
}

void	Orifice::fnSetState(double mdotv)
{
	VState(1,mdotv);
}

bool	Orifice::fnReEvalState(double& mdotv)
{
	char errstr[512];
	//Capture boundary conditions
	double P1, P2, h3;
	P1 = VINPUTS(1);
	P2 = VINPUTS(2);
	h3 = VINPUTS(3);

	//Capture geometry
	double A, C_d;
	A	= VGeometry(1);
	C_d	= VGeometry(2);


	//Re-calculate state variable
	double dP, v3;
	r134astate	inlet;
	if(!inlet.setstate(P2,h3))
	{
		sprintf(errstr,"Could not set refrigerant properties with %fkPa and %fkJ/kg",P2,h3);
		errorlog.Add("Orifice::fnReEvalState ",errstr);
		return false;
	}
	v3	=	inlet.getv();
	dP		=	(P2 - P1)*1e3;	//convert to Pa, for units consistency
	if(dP>=0)
		mdotv	=	C_d*A*sqrt(dP/v3);
	else
		mdotv	=	0.0;

	return true;
}

bool	Orifice::fnAdvance1sec()
{
	double mdotcl;
	if(!fnReEvalState(mdotcl))	return false;
	VState(1,mdotcl);
	VOUTPUTS(1,mdotcl);
	return true;
}

//*******************************************************************************************

VapCompCentLiqChiller::VapCompCentLiqChiller()
{

	iN_INPUTS	=	5;
	iN_OUTPUTS	=	28;

	VINPUTS.make(iN_INPUTS);
	VOUTPUTS.make(iN_OUTPUTS);

	dTIME	=	0.0;

	//Load property tables
	PROP_TABLES_LOADED = loadtables(R134A,1);
	if(!PROP_TABLES_LOADED)
		errorlog.Add("All or some of the property tables could not be loaded");

	VState.make(6);								///////////////////////CHANGED//////////////////////////
	Evaporator	= new FVShellTubeHX;
	Condenser	= new FVShellTubeHX;
	Compressor	= new CentComp;
	Valve		= new ThermoStaticExValve;
	CoolingLineBypass		= new Orifice;
	bINITIALIZED = false;
	bNONCONDENSABLES = false;
}

bool	VapCompCentLiqChiller::fnInitialized()
{
	return bINITIALIZED;
}

void	VapCompCentLiqChiller::fnSetNCMass(double NCMass)
{
	dNCMass = NCMass;
	bNONCONDENSABLES = true;
}

double	VapCompCentLiqChiller::fnGetNCMass()
{
	return dNCMass;
}

bool	VapCompCentLiqChiller::fnDefineGeometry()
{
	char		geomfile[]	=	"Geometries\\SystemGeometry.txt";
	ifstream	gfile(geomfile);
	if(!gfile)	return false;

	char *evapgeomfile, *condgeomfile, *valvgeomfile, *compgeomfile, *orifgeomfile;
	char line[512];

	gfile >> line;
	evapgeomfile = new char[strlen(line)];
	strcpy(evapgeomfile,line);

	gfile >> line;
	condgeomfile = new char[strlen(line)];
	strcpy(condgeomfile,line);

	gfile >> line;
	compgeomfile = new char[strlen(line)];
	strcpy(compgeomfile,line);

	gfile >> line;
	valvgeomfile = new char[strlen(line)];
	strcpy(valvgeomfile,line);

	gfile >> line;
	orifgeomfile = new char[strlen(line)];
	strcpy(orifgeomfile,line);

	if(!Evaporator->fnDefineGeometry(evapgeomfile))
	{
		errorlog.Add("VapCompCentLiqChiller::fnDefineGeometry ","Evaporator construction failed.");
		return false;
	}
	if(!Condenser->fnDefineGeometry(condgeomfile))
	{
		errorlog.Add("VapCompCentLiqChiller::fnDefineGeometry ","Condenser construction failed.");
		return false;
	}
	if(!Compressor->fnDefineGeometry(compgeomfile))
	{
		errorlog.Add("VapCompCentLiqChiller::fnDefineGeometry ","Compressor construction failed.");
		return false;
	}
	if(!Valve->fnDefineGeometry(valvgeomfile))
	{
		errorlog.Add("VapCompCentLiqChiller::fnDefineGeometry ","Valve construction failed.");
		return false;
	}
	if(!CoolingLineBypass->fnDefineGeometry(orifgeomfile))
	{
		errorlog.Add("VapCompCentLiqChiller::fnDefineGeometry ","Orifice construction failed.");
		return false;
	}

	return true;
}

bool	VapCompCentLiqChiller::fnDefineGeometry(const char* sysgeomfile)
{
	ifstream	gfile(sysgeomfile);
	if(!gfile)	return false;

	char *evapgeomfile, *condgeomfile, *valvgeomfile, *compgeomfile, *orifgeomfile;
	char line[512];

	gfile >> line;
	evapgeomfile = new char[strlen(line)];
	strcpy(evapgeomfile,line);

	gfile >> line;
	condgeomfile = new char[strlen(line)];
	strcpy(condgeomfile,line);

	gfile >> line;
	compgeomfile = new char[strlen(line)];
	strcpy(compgeomfile,line);

	gfile >> line;
	valvgeomfile = new char[strlen(line)];
	strcpy(valvgeomfile,line);

	gfile >> line;
	orifgeomfile = new char[strlen(line)];
	strcpy(orifgeomfile,line);

	if(!Evaporator->fnDefineGeometry(evapgeomfile))
	{
		errorlog.Add("VapCompCentLiqChiller::fnDefineGeometry ","Evaporator construction failed.");
		return false;
	}
	if(!Condenser->fnDefineGeometry(condgeomfile))
	{
		errorlog.Add("VapCompCentLiqChiller::fnDefineGeometry ","Condenser construction failed.");
		return false;
	}
	if(!Compressor->fnDefineGeometry(compgeomfile))
	{
		errorlog.Add("VapCompCentLiqChiller::fnDefineGeometry ","Compressor construction failed.");
		return false;
	}
	if(!Valve->fnDefineGeometry(valvgeomfile))
	{
		errorlog.Add("VapCompCentLiqChiller::fnDefineGeometry ","Valve construction failed.");
		return false;
	}
	if(!CoolingLineBypass->fnDefineGeometry(orifgeomfile))
	{
		errorlog.Add("VapCompCentLiqChiller::fnDefineGeometry ","Orifice construction failed.");
		return false;
	}

	return true;
}

bool	VapCompCentLiqChiller::fnSetState(double time, double evapP, Matrix& evapS, 
							 double condP, Matrix& condS, 
							 int mode, Vector& compS, 
							 double mdotv, double Tbulb, double mdotcl,double Tewo,double NCMass)
{
	if(!Evaporator->fnSetState(time,evapP,evapS))
	{
		errorlog.Add("VapCompCentLiqChiller::fnSetState ","Evaporator initialization failed." );
		return false;
	}
	if(!Condenser->fnSetState(time,condP,condS))
	{
		errorlog.Add("VapCompCentLiqChiller::fnSetState ","Condenser initialization failed." );
		return false;
	}
	if(!Compressor->fnSetState(time,mode,compS))
	{
		errorlog.Add("VapCompCentLiqChiller::fnSetState ","Compressor initialization failed." );
		return false;
	}
	if(!Valve->fnSetState(time,mdotv,Tbulb))		
	{
		errorlog.Add("VapCompCentLiqChiller::fnSetState ","Valve initialization failed." );
		return false;
	}
	CoolingLineBypass->fnSetState(mdotcl);

	VState(1,evapP);	VState(2,condP);
	VState(3,evapS(evapS.rows(),1));
	VState(4,condS(condS.rows(),1));
	VState(5,Tewo);

	bINITIALIZED = true;
	if(NCMass>0.0)
	{
		fnSetNCMass(NCMass);
		bNONCONDENSABLES = true;
	}

	return true;
}


bool	VapCompCentLiqChiller::fnSetState(double Tewo, double Tcwo, double Mtotal, double NCMass)
{
	char errstr[512];
	r134astate	state;
	double time = 0.0;
	double Pevap, Mevap;
	double Pcond, hcond, vcond, Mcond, Vcond;

	if(!state.setstateT(Tewo))
	{
		sprintf(errstr,"Could not set refrigerant properties with %fC",Tewo);
		errorlog.Add("VapCompCentLiqChiller::fnSetState ",errstr);
		return false;
	}

	Pevap	=	state.getP();
	Pcond	=	Pevap;
	if(!state.setstatePT(Pcond,Tcwo))
	{
		sprintf(errstr,"Could not set refrigerant properties with %fkPa and %fC",Pcond,Tcwo);
		errorlog.Add("VapCompCentLiqChiller::fnSetState ",errstr);
		return false;
	}
	vcond	=	state.getv();
	hcond	=	state.geth();
	Vcond	=	Condenser->Volume();
	Mcond	=	Vcond/vcond;
	Mevap	=	Mtotal - Mcond;

	if(!Evaporator->fnSetState(time,Pevap,Mevap))
	{
		errorlog.Add("VapCompCentLiqChiller::fnSetState ","Evaporator initialization failed." );
		return false;
	}
	if(!Condenser->fnSetState(time,Pcond,Mcond))
	{
		errorlog.Add("VapCompCentLiqChiller::fnSetState ","Condenser initialization failed." );
		return false;
	}

	if(!Compressor->fnSetState(hcond))
	{
		errorlog.Add("VapCompCentLiqChiller::fnSetState ","Compressor initialization failed." );
		return false;
	}
	if(!Valve->fnSetState(Tewo))
	{
		errorlog.Add("VapCompCentLiqChiller::fnSetState ","Valve initialization failed." );
		return false;
	}
	CoolingLineBypass->fnSetState(0.0);

	double h1, h3;
	h1 = state.gethl();
	h3 = state.gethv();
	VState(1,Pevap);	VState(2,Pcond);	VState(3,h1);	VState(4,h3);	VState(5,Tewo);
	VState(6,h1);			///////////////////////ADDED//////////////////////////

	bINITIALIZED = true;
	if(NCMass>0.0)
	{
		fnSetNCMass(NCMass);
		bNONCONDENSABLES = true;
	}

	return true;
}

bool	VapCompCentLiqChiller::fnSaveState()
{
	ofstream	sfile("IOFiles\\SystemState.txt");
	if(!sfile)	return false;

	double	Pevap;
	Matrix	EvapS;
	Evaporator->fnGetState(Pevap,EvapS);

	double	Pcond;
	Matrix	CondS;
	Condenser->fnGetState(Pcond,CondS);

	int	mode;
	Vector CompS;
	Compressor->fnGetState(mode,CompS);

	Vector ValvS;
	Valve->fnGetState(ValvS);

	int enodes, cnodes;
	enodes = EvapS.rows();
	cnodes = CondS.rows();

	sfile << Condenser->dTIME << endl;

	sfile << Pevap << "\t" << Pcond << "\t" << EvapS(enodes,1)	<< "\t" << CondS(cnodes,1) << "\t" << VState(5) << endl;
	sfile << endl;
	for(int i=1;i<=enodes;i++)
		//Evap....hi.....................Tti...................Twi..................mri..................Qr
		sfile << EvapS(i,1) << "\t" << EvapS(i,2) << "\t" << EvapS(i,3) << "\t" << EvapS(i,4) << "\t" << EvapS(i,5) << endl;
	sfile << endl;
	for(int i=1;i<=cnodes;i++)
		//Cond....hi.....................Tti...................Twi..................mri..................Qr
		sfile << CondS(i,1) << "\t" << CondS(i,2) << "\t" << CondS(i,3) << "\t" << CondS(i,4) << "\t" << CondS(i,5) << endl;
	sfile << endl;
	//Comp...mode....h2.................mdotc...............Power................Qloss...............gamma
	sfile << mode << "\t" << CompS(1) << "\t" << CompS(2) << "\t" << CompS(3) << "\t" << CompS(4) << "\t" << CompS(5) << endl;
	sfile << endl;
	//Valve..Tbulb.......mdotv
	sfile << ValvS(2) <<  "\t" << ValvS(1) << endl;

	sfile << fnGetNCMass() << endl;

	sfile.close();

	return true;
}

bool	VapCompCentLiqChiller::fnSaveState(const char* statefile)
{
	ofstream	sfile(statefile);
	if(!sfile)	return false;

	double	Pevap;
	Matrix	EvapS;
	Evaporator->fnGetState(Pevap,EvapS);

	double	Pcond;
	Matrix	CondS;
	Condenser->fnGetState(Pcond,CondS);

	int	mode;
	Vector CompS;
	Compressor->fnGetState(mode,CompS);

	Vector ValvS;
	Valve->fnGetState(ValvS);

	int enodes, cnodes;
	enodes = EvapS.rows();
	cnodes = CondS.rows();

	sfile << Condenser->dTIME << endl;

	sfile << Pevap << "\t" << Pcond << "\t" << EvapS(enodes,1)	<< "\t" << CondS(cnodes,1) << "\t" << VState(5) << endl;
	sfile << endl;
	for(int i=1;i<=enodes;i++)
		//Evap....hi.....................Tti...................Twi..................mri..................Qr
		sfile << EvapS(i,1) << "\t" << EvapS(i,2) << "\t" << EvapS(i,3) << "\t" << EvapS(i,4) << "\t" << EvapS(i,5) << endl;
	sfile << endl;
	for(int i=1;i<=cnodes;i++)
		//Cond....hi.....................Tti...................Twi..................mri..................Qr
		sfile << CondS(i,1) << "\t" << CondS(i,2) << "\t" << CondS(i,3) << "\t" << CondS(i,4) << "\t" << CondS(i,5) << endl;
	sfile << endl;
	//Comp...mode....h2.................mdotc...............Power................Qloss...............gamma
	sfile << mode << "\t" << CompS(1) << "\t" << CompS(2) << "\t" << CompS(3) << "\t" << CompS(4) << "\t" << CompS(5) << endl;
	sfile << endl;
	//Valve..Tbulb.......mdotv
	sfile << ValvS(2) <<  "\t" << ValvS(1) << endl;

	sfile << fnGetNCMass() << endl;

	sfile.close();

	return true;
}

bool	VapCompCentLiqChiller::fnLoadState()
{
	char errstr[512];
	ifstream	ifile("IOFiles\\Initial_FULL.txt");
	if(!ifile)	return false;

	double Pevap, Pcond, h1, h3;

	Matrix EvapS;
	int enodes;
	enodes	=	Evaporator->NODES();
	EvapS.make(enodes,5);
	double Tt, Tw, Q, mr;

	Matrix CondS;
	int cnodes;
	cnodes	=	Condenser->NODES();
	CondS.make(cnodes,5);

	Vector	CompS;
	CompS.make(5);
	int mode;
	double	h2, mdotc, Pwr, Qloss, gamma;

	double Tbulb;


	dTIME = 0.0;
	ifile >> Pevap >> h1;
	r134astate state;
	if(!state.setstate(Pevap,h1))
	{
		sprintf(errstr,"Could not set refrigerant properties with %fkPa and %fkJ/kg",Pevap,h1);
		errorlog.Add("VapCompCentLiqChiller::fnLoadState ",errstr);
		return false;
	}
	Tt = state.getT();
	Tw = state.getT();
	Tbulb = Tt;
	mr = 0.0;
	Q = 0.0;
	for(int i=1;i<=enodes;i++)
	{
		EvapS(i,1,h1);	EvapS(i,2,Tt);	EvapS(i,3,Tw);	EvapS(i,4,mr);	EvapS(i,5,Q);
	}

	ifile >> Pcond >> h3;
	if(!state.setstate(Pcond,h3))
	{
		sprintf(errstr,"Could not set properties with %fkPa and %fkJ/kg",Pcond,h3);
		errorlog.Add("VapCompCentLiqChiller::fnLoadState ",errstr);
		return false;
	}
	Tt = state.getT();
	Tw = state.getT();
	for(int i=1;i<=cnodes;i++)
	{
		CondS(i,1,h3);	CondS(i,2,Tt);	CondS(i,3,Tw);	CondS(i,4,mr);	CondS(i,5,Q);
	}

	mode = STARTUP;
	h2 = h3;
	mdotc = 0.0;
	Pwr = 0.0;
	Qloss = 0.0;
	gamma = 0.05;
	CompS(1,h2);	CompS(2,mdotc);	CompS(3,Pwr);	CompS(4,Qloss);	CompS(5,gamma);

	if(!Evaporator->fnSetState(dTIME,Pevap,EvapS))
	{
		errorlog.Add("VapCompCentLiqChiller::fnLoadState ","Evaporator initialization failed");
		return false;
	}

	if(!Condenser->fnSetState(dTIME,Pcond,CondS))
	{
		errorlog.Add("VapCompCentLiqChiller::fnLoadState ","Condenser initialization failed");
		return false;
	}

	if(!Compressor->fnSetState(dTIME,mode,CompS))
	{
		errorlog.Add("VapCompCentLiqChiller::fnLoadState ","Compressor initialization failed");
		return false;
	}
	if(!Valve->fnSetState(dTIME,0.0,Tbulb))
	{
		errorlog.Add("VapCompCentLiqChiller::fnLoadState ","Valve initialization failed");
		return false;
	}

	VState(1,Pevap);	VState(2,Pcond);	VState(3,h1);	VState(4,h3);	VState(5,Tw);
	VState(6,EvapS(BULBLOCN,1));		///////////////////////ADDED//////////////////////////
	bINITIALIZED = true;

	double NCMass;
	ifile >> NCMass;
	if(NCMass>0.0)
	{
		fnSetNCMass(NCMass);
		bNONCONDENSABLES = true;
	}

	ifile.close();
	return true;
}

bool	VapCompCentLiqChiller::fnLoadState(const char* sysstate)
{
	char errstr[512];
	ifstream	ifile(sysstate);
	if(!ifile)	return false;

	char line[512];
	double Pevap, Pcond, h1, h3, Tewo;

	Matrix EvapS;
	int enodes;
	enodes	=	Evaporator->NODES();
	EvapS.make(enodes,5);
	double h, Tt, Tw, Q, mr;

	Matrix CondS;
	int cnodes;
	cnodes	=	Condenser->NODES();
	CondS.make(cnodes,5);

	Vector	CompS;
	CompS.make(5);
	int mode;
	double	h2, mdotc, Pwr, Qloss, gamma;

	double Tbulb;

	if(!strcmp(sysstate,"IOFiles\\Initial_MINIMAL.txt"))
		return fnLoadState(MINIMAL);
	else if(!strcmp(sysstate,"IOFiles\\Initial_FULL.txt"))
	{
		dTIME = 0.0;
		ifile.getline(line,512,'\n');
		if(sscanf(line,"%lf%lf",&Pevap,&h1)!=2)
		{
			errorlog.Add("VapCompCentLiqChiller::fnLoadState ","Incorrect file format.");
			return false;
		}
		r134astate state;
		if(!state.setstate(Pevap,h1))
		{
			sprintf(errstr,"Could not set refrigerant properties with %fkPa and %fkJ/kg",Pevap,h1);
			errorlog.Add("VapCompCentLiqChiller::fnLoadState ",errstr);
			return false;
		}
		Tt = state.getT();
		Tw = state.getT();
		Tbulb = Tt;
		mr = 0.0;
		Q = 0.0;
		for(int i=1;i<=enodes;i++)
		{
			EvapS(i,1,h1);	EvapS(i,2,Tt);	EvapS(i,3,Tw);	EvapS(i,4,mr);	EvapS(i,5,Q);
		}

		ifile.getline(line,512,'\n');
		if(sscanf(line,"%lf%lf",&Pcond,&h3)!=2)
		{
			errorlog.Add("VapCompCentLiqChiller::fnLoadState ","Incorrect file format.");
			return false;
		}
		if(!state.setstate(Pcond,h3))
		{
			sprintf(errstr,"Could not set refrigerant properties with %fkPa and %fkJ/kg",Pcond,h3);
			errorlog.Add("VapCompCentLiqChiller::fnLoadState ",errstr);
			return false;
		}
		Tt = state.getT();
		Tw = state.getT();
		for(int i=1;i<=cnodes;i++)
		{
			CondS(i,1,h3);	CondS(i,2,Tt);	CondS(i,3,Tw);	CondS(i,4,mr);	CondS(i,5,Q);
		}

		mode = STARTUP;
		h2 = h3;
		mdotc = 0.0;
		Pwr = 0.0;
		Qloss = 0.0;
		gamma = 0.05;
		CompS(1,h2);	CompS(2,mdotc);	CompS(3,Pwr);	CompS(4,Qloss);	CompS(5,gamma);

		double NCMass;
		ifile >> NCMass;
		if(NCMass>0.0)
		{
			fnSetNCMass(NCMass);
			bNONCONDENSABLES = true;
		}
	}
	else
	//i.e. Input file is NOT Initial_FULL.txt or Initial_MINIMAL.txt
	{
		ifile.getline(line,512,'\n');
		if(sscanf(line,"%lf",&dTIME)!=1)
		{
			errorlog.Add("VapCompCentLiqChiller::fnLoadState ","Incorrect file format.");
			return false;
		}
		ifile.getline(line,512,'\n');
		if(sscanf(line,"%lf%lf%lf%lf%lf",&Pevap,&Pcond,&h1,&h3,&Tewo)!=5)
		{
			errorlog.Add("VapCompCentLiqChiller::fnLoadState ","Incorrect file format.");
			return false;
		}
		ifile.getline(line,512,'\n');
		VState(1,Pevap);	VState(2,Pcond);	VState(3,h1);	VState(4,h3);	VState(5,Tewo);

		for(int i=1;i<=enodes;i++)
		{
			ifile.getline(line,512,'\n');
			if(sscanf(line,"%lf%lf%lf%lf%lf",&h,&Tt,&Tw,&mr,&Q)!=5)
			{
				errorlog.Add("VapCompCentLiqChiller::fnLoadState ","Incorrect file format.");
				return false;
			}
			EvapS(i,1,h);	EvapS(i,2,Tt);	EvapS(i,3,Tw);	EvapS(i,4,mr);	EvapS(i,5,Q);
		}
		ifile.getline(line,512,'\n');

		for(int i=1;i<=cnodes;i++)
		{
			ifile.getline(line,512,'\n');
			if(sscanf(line,"%lf%lf%lf%lf%lf",&h,&Tt,&Tw,&mr,&Q)!=5)
			{
				errorlog.Add("VapCompCentLiqChiller::fnLoadState ","Incorrect file format.");
				return false;
			}
			CondS(i,1,h);	CondS(i,2,Tt);	CondS(i,3,Tw);	CondS(i,4,mr);	CondS(i,5,Q);
		}
		ifile.getline(line,512,'\n');

		ifile.getline(line,512,'\n');
		if(sscanf(line,"%d%lf%lf%lf%lf%lf",&mode,&h2,&mdotc,&Pwr,&Qloss,&gamma)!=6)
		{
			errorlog.Add("VapCompCentLiqChiller::fnLoadState ","Incorrect file format.");
			return false;
		}
		CompS(1,h2);	CompS(2,mdotc);	CompS(3,Pwr);	CompS(4,Qloss);	CompS(5,gamma);
		ifile.getline(line,512,'\n');
		ifile.getline(line,512,'\n');
		if(sscanf(line,"%lf",&Tbulb)!=1)
		{
			errorlog.Add("VapCompCentLiqChiller::fnLoadState ","Incorrect file format.");
			return false;
		}
		ifile.getline(line,512,'\n');
		ifile.getline(line,512,'\n');
		double NCMass;
		if(sscanf(line,"%lf",&NCMass)!=1)
		{
			errorlog.Add("VapCompCentLiqChiller::fnLoadState ","Incorrect file format.");
			return false;
		}
		if(NCMass>0.0)
		{
			fnSetNCMass(NCMass);
			bNONCONDENSABLES = true;
		}
	}


	if(!Evaporator->fnSetState(dTIME,Pevap,EvapS))
	{
		errorlog.Add("VapCompCentLiqChiller::fnLoadState ","Evaporator initialization failed");
		return false;
	}

	if(!Condenser->fnSetState(dTIME,Pcond,CondS))
	{
		errorlog.Add("VapCompCentLiqChiller::fnLoadState ","Condenser initialization failed");
		return false;
	}

	if(!Compressor->fnSetState(dTIME,mode,CompS))
	{
		errorlog.Add("VapCompCentLiqChiller::fnLoadState ","Compressor initialization failed");
		return false;
	}
	if(!Valve->fnSetState(dTIME,0.0,Tbulb))
	{
		errorlog.Add("VapCompCentLiqChiller::fnLoadState ","Valve initialization failed");
		return false;
	}

	VState(1,Pevap);	VState(2,Pcond);	VState(3,h1);	VState(4,h3);	VState(5,Tw);
	VState(6,EvapS(BULBLOCN,1));		///////////////////////ADDED//////////////////////////

	bINITIALIZED = true;

	ifile.close();
	return true;
}

bool	VapCompCentLiqChiller::fnLoadState(int Choice)
{
	if(Choice!=MINIMAL&&Choice!=FULL)
	{
		errorlog.Add("VapCompCentLiqChiller::fnLoadState ","Invalid Choice to fnLoadState.",(double)Choice);
		return false;
	}
	if(Choice==MINIMAL)
	{
		char	sysstate[]	=	"IOFiles\\Initial_MINIMAL.txt";
		ifstream ifile(sysstate);
		if(!ifile)
		{
			errorlog.Add("VapCompCentLiqChiller::fnLoadState ","File not found.");
			return false;
		}
		double Tewo, Tcwo, Mtotal;
		char line[512]; int count;
		ifile.getline(line,512,'\n');
		count = sscanf(line,"%lf%lf%lf",&Tewo,&Tcwo,&Mtotal);
		if(count!=3)
		{
			errorlog.Add("VapCompCentLiqChiller::fnLoadState ","File format is not correct.");
			return false;
		}
		double NCMass;
		ifile >> NCMass;
		if(!fnSetState(Tewo,Tcwo,Mtotal,NCMass))
		{
			errorlog.Add("VapCompCentLiqChiller::fnLoadState ","Could not do minimal initialization.");
			return false;
		}
		ifile.close();
	}
	if(Choice==FULL)
	{
		if(!fnLoadState("IOFiles\\Initial_FULL.txt"))
		{
			errorlog.Add("VapCompCentLiqChiller::fnLoadState ","Could not do full initialization.");
			return false;
		}
	}

	return true;
}

bool	VapCompCentLiqChiller::fnAdvance1sec(int METHOD)
{
	char errstr[512];

	//Capture current state;
	double Pevap, Pcond, h1, h3, Tewo, hbulb;
	Pevap	=	VState(1);	Pcond	=	VState(2);
	h1		=	VState(3);	h3		=	VState(4);
	Tewo	=	VState(5);	hbulb	=	VState(6);

	//Capture system bcs
	double Tewin, Tcwin, TewoSet, mewat, mcwat;
	Tewin	=	VINPUTS(1);	Tcwin	=	VINPUTS(2);	TewoSet	=	VINPUTS(3);
	mewat	=	VINPUTS(4);	mcwat	=	VINPUTS(5);

	//COMPRESSOR
	Compressor->VINPUTS(1,Pevap);	Compressor->VINPUTS(2,h1);	Compressor->VINPUTS(3,Pcond);
	Compressor->VINPUTS(4,TewoSet);	Compressor->VINPUTS(5,Tewo);
	if(!Compressor->fnAdvance1sec())
	{
		sprintf(errstr,"Compressor execution failed at time %f",dTIME);
		errorlog.Add("VapCompCentLiqChiller::fnAdvance1sec ",errstr);
		return false;
	}
	double h2, mdotc, Pwr, Qloss,cmpengyIMB;
	h2		=	Compressor->VOUTPUTS(1);
	mdotc	=	Compressor->VOUTPUTS(2);
	Pwr		=	Compressor->VOUTPUTS(3);
	Qloss	=	Compressor->VOUTPUTS(4);
	cmpengyIMB	=	Pwr - Qloss - mdotc*(h2-Compressor->VINPUTS(2));

	//VALVE
	Valve->VINPUTS(1,Pevap);	Valve->VINPUTS(2,hbulb);	Valve->VINPUTS(3,Pcond);	Valve->VINPUTS(4,h3);
	if(!Valve->fnAdvance1sec(Tewo))
	{
		sprintf(errstr,"Valve execution failed at time %f",dTIME);
		errorlog.Add("VapCompCentLiqChiller::fnAdvance1sec ",errstr);
		return false;
	}
	double mdotv, lift, Area;
	mdotv	=	Valve->VOUTPUTS(1);
	lift	=	Valve->VOUTPUTS(4);
	Area	=	Valve->VOUTPUTS(5);

	//ORIFICE
	double mdotcl;	double mvalve;
	CoolingLineBypass->VINPUTS(1,Pevap);	CoolingLineBypass->VINPUTS(2,Pcond);	CoolingLineBypass->VINPUTS(3,h3);
	if(!CoolingLineBypass->fnAdvance1sec())
	{
		sprintf(errstr,"Cooling line orifice execution failed at time %f",dTIME);
		errorlog.Add("VapCompCentLiqChiller::fnAdvance1sec ",errstr);
		return false;
	}
	mdotcl	=	CoolingLineBypass->VOUTPUTS(1);
	mvalve	=	mdotv+mdotcl;

	//CONDENSER
	Condenser->VINPUTS(1,mcwat);	Condenser->VINPUTS(2,Tcwin);	Condenser->VINPUTS(3,mdotc);
	Condenser->VINPUTS(4,mvalve);	Condenser->VINPUTS(5,h2);
	if(!Condenser->fnAdvance1sec(METHOD))
	{
		sprintf(errstr,"Condenser execution failed at time %f",dTIME);
		errorlog.Add("VapCompCentLiqChiller::fnAdvance1sec ",errstr);
		return false;
	}
	double Tcwo, Qc, Tsub, cmassIMB, h4,CREFQTY,cengyIMB;
	Pcond	=	Condenser->VOUTPUTS(1);
	h3		=	Condenser->VOUTPUTS(2);
	Tcwo	=	Condenser->VOUTPUTS(3);
	Qc		=	Condenser->VOUTPUTS(4);
	Tsub	=	Condenser->VOUTPUTS(5);
	cengyIMB=	Condenser->VOUTPUTS(6);
	CREFQTY	=	Condenser->VOUTPUTS(7);
	cmassIMB=	Condenser->VOUTPUTS(8);
	h4		=	(mdotv+mdotcl<MYZERO)?h3:(h3+Qloss/(mdotv+mdotcl));

	//Non Condensables correction
	if(bNONCONDENSABLES)
	{
		double Pnc, Vhx, Mnc, Rnc, Thx;
		r134astate ncstate;
		ncstate.setstateP(Pcond);
		Thx = ncstate.getTs() + 273.16;		//K
		Rnc = 0.2968;	//kJ/kg-K
		Mnc = fnGetNCMass();
		Vhx = Condenser->Volume();

		Pnc = (Rnc*Thx*Mnc)/Vhx;

		Pcond = Pcond + Pnc;

	}

	//EVAPORATOR
	Evaporator->VINPUTS(1,mewat);	Evaporator->VINPUTS(2,Tewin);	Evaporator->VINPUTS(3,mvalve);
	Evaporator->VINPUTS(4,mdotc);	Evaporator->VINPUTS(5,h4);
	if(!Evaporator->fnAdvance1sec(METHOD))
	{
		sprintf(errstr,"Evaporator execution failed at time %f",dTIME);
		errorlog.Add("VapCompCentLiqChiller::fnAdvance1sec ",errstr);
		return false;
	}
	double Qe, Tsh, emassIMB,EREFQTY,eengyIMB;
	Pevap	=	Evaporator->VOUTPUTS(1);
	h1		=	Evaporator->VOUTPUTS(2);
	Tewo	=	Evaporator->VOUTPUTS(3);
	Qe		=	Evaporator->VOUTPUTS(4);
	Tsh		=	Evaporator->VOUTPUTS(5);
	eengyIMB=	Evaporator->VOUTPUTS(6);
	EREFQTY	=	Evaporator->VOUTPUTS(7);
	emassIMB=	Evaporator->VOUTPUTS(8);
	hbulb	=	Evaporator->VOUTPUTS(9);

	VState(1,Pevap);	VState(2,Pcond);	VState(3,h1);
	VState(4,h3);		VState(5,Tewo);		VState(6,hbulb);

	dTIME	+=	1.0;

	VOUTPUTS(1,Pevap);		VOUTPUTS(2,Pcond);		
	VOUTPUTS(3,mdotc);		VOUTPUTS(4,mdotv);		VOUTPUTS(5,mdotcl);		VOUTPUTS(6,mvalve);		
	VOUTPUTS(7,Pwr);		VOUTPUTS(8,Qloss);		VOUTPUTS(9,Qc);			VOUTPUTS(10,Qe);
	VOUTPUTS(11,Tewo);		VOUTPUTS(12,Tcwo);		VOUTPUTS(13,Tsh);		VOUTPUTS(14,Tsub);
	VOUTPUTS(15,cmassIMB);	VOUTPUTS(16,emassIMB);
	VOUTPUTS(17,cmpengyIMB);VOUTPUTS(18,cengyIMB);	VOUTPUTS(19,eengyIMB);
	VOUTPUTS(20,h1);		VOUTPUTS(21,h2);		VOUTPUTS(22,h3);		VOUTPUTS(23,h4);
	VOUTPUTS(24,lift);		VOUTPUTS(25,Area);
	VOUTPUTS(26,CREFQTY);	VOUTPUTS(27,EREFQTY);	VOUTPUTS(28,CREFQTY+EREFQTY);

	return true;
}

void	VapCompCentLiqChiller::shutoff()
{
	Compressor->shutoff();
}

/******************************************************************************************/
/****************M O V I N G * B O U N D A R Y * F O R M U L A T I O N*********************/
/****************************************O F***********************************************/
/***************S H E L L * - * T U B E * H E A T * E X C H A N G E R S********************/
/******************************************************************************************/

/************************************************************************/
/**********************C O N D E N S E R*********************************/
/************************************************************************/

MBShellTubeCondenser::MBShellTubeCondenser()
{
	VINPUTS.make(5);
	VOUTPUTS.make(12);
	VGeometry.make(10);
	VState.make(12);
	rpropertylist.make(20);
	wpropertylist.make(5);
}//MBShellTubeCondenser::MBShellTubeCondenser

bool	MBShellTubeCondenser::fnDefineGeometry(const char* fname)
{
	ifstream	infile;
	infile.open(fname);

	int N, scrap;
	double tubeid, tubeod, tubelen, tubeCp, tuberho, foulingfactor, shellid;

	infile >> iHXTYPE >> scrap >> N;
	infile >> tubeid >> tubeod >> tubelen >> tubeCp >> tuberho >> foulingfactor >> shellid;
	infile.close();

	NTUBES = N;
	VGeometry(1,tubeid);	VGeometry(2,tubeod);	VGeometry(3,tubelen);

	double TotShellVolume, TotHTArea, TotTubeCap, TotWatCap, farea;
	double waterrho = 995.0, waterCp = 4.1868;

	farea			=	(MYPI/4)*(shellid*shellid-NTUBES*tubeod*tubeod);
	TotShellVolume	=	TotVolume = farea;		//Volume per m
	TotHTArea		=	NTUBES*MYPI*tubeod;		//Outside ht area per m
	TotTubeCap		=	NTUBES*(MYPI/4)*(tubeod*tubeod-tubeid*tubeid)*tuberho*tubeCp;	//Tube capacitance per m
	TotWatCap		=	NTUBES*(MYPI/4)*tubeid*tubeid*waterrho*waterCp;					//Water capacitance per m

	VGeometry(4,farea);
	VGeometry(5,TotHTArea);
	VGeometry(6,TotHTArea*tubeid/tubeod);
	VGeometry(7,TotTubeCap);
	VGeometry(8,double(NTUBES));
	VGeometry(9,TotWatCap);
	VGeometry(10,foulingfactor);

	return true;
}//MBShellTubeCondenser::fnDefineGeometry

bool	MBShellTubeCondenser::fnInitialize(Vector& init)
{
	iN_INPUTS = 5;	iN_OUTPUTS = 12;

	VState(1,init(1));		//Pressure, kPa
	VState(2,init(2));		//L1, m
	VState(3,init(3));		//L2, m
	VState(4,init(4));		//hout, m
	VState(5,init(5));		//mdotA, kg/s
	VState(6,init(6));		//mdotB, kg/s	
	VState(7,init(7));		//Tt1, C
	VState(8,init(8));		//Tt2, C
	VState(9,init(9));		//Tt3, C
	VState(10,init(10));	//Tw1, C
	VState(11,init(11));	//Tw2, C
	VState(12,init(12));	//Tw3, C

	double	L1, L2, L;
	L1	=	init(2);
	L2	=	init(3);
	L = VGeometry(3);

	if(L1 < 0.0 || L1 < 0.0 || L1 > L || L2 > L)
	{
		errorlog.Add("MBShelltubeCondenser::fnInitialize: ","Invalid initialization.");
		return false;
	}

	r134astate	initstate;
	double P, hout, hv, hl;
	P = init(1);
	hout = init(4);
	if(!initstate.setstateP(P))
	{
		errorlog.Add("MBShellTubecondenser::fnInitialize: ","Could  not set initial state.");
		return false;
	}

	hv = initstate.gethv();	hl = initstate.gethl();

	if(hout > hv)
		iOPMODE	=	SH;
	else if(hout < hl)
		iOPMODE	=	SHTPSC;
	else
		iOPMODE	=	SHTP;

	double rho, farea, u;
	farea = VGeometry(4);

	initstate.setstate(P,hout);
	rho	=	initstate.getrho();
	u	=	initstate.getu();

	dORIGREFQTY =	rho*L*farea;
	dORIGREFENGY=	rho*u*L*farea;

	dTIME	=	0.0;
	dTSTEP	=	0.01;
	dTOTMASSIN	=	0.0;
	dTOTMASSOUT	=	0.0;
	dTOTENGYIN	=	0.0;
	dTOTENGYOUT	=	0,0;

	dHRINOLD.make(HINDEPTH);

	return true;

}//MBShellTubeCondenser::fnInitialize

double	MBShellTubeCondenser::whtcoeff(double Tw, double vel)
{
	watstate water;
	water.setstate(Tw);

	wpropertylist(1,water.getmu());
	wpropertylist(2,water.getk());
	wpropertylist(3,water.getCp());
	wpropertylist(4,water.getrho());
	wpropertylist(5,vel);

	return fnWatHTCoeff(wpropertylist);
}//MBShellTubeCondenser::whtcoeff

double	MBShellTubeCondenser::fnWatHTCoeff(Vector& wpropertylist)
{

	double	Pr,Re,Nu,rho,mu,Cp,k,h,f,vel,tubeid;
	mu	=	wpropertylist(1);
	k	=	wpropertylist(2);
	Cp	=	wpropertylist(3);
	rho	=	wpropertylist(4);
	vel	=	wpropertylist(5);

	tubeid = VGeometry(1);

	Pr	= mu*Cp*1e3/k;
	Re	= rho*vel*tubeid/mu;

	f	=	1/pow(0.79*log(Re)-1.64,2);
	Nu	=	(f*Re*Pr/8)/(1.07+12.7*sqrt(f/8)*(pow(Pr,0.67)-1));
	h	=	2.4*k*1e-3*Nu/tubeid;

	h	=	3.0*h;
	return h;
}//MBShellTubeCondenser::fnWatHTCoeff

double	MBShellTubeCondenser::rhtcoeff(double P, double h, double vel, double Tt)
{
	r134astate state;
	state.setstate(P,h);

	rpropertylist(1,state.getTs());		rpropertylist(3,state.getmu());		rpropertylist(4,state.getCp());
	rpropertylist(5,state.getk());		rpropertylist(6,state.getrho());	rpropertylist(7,state.getphase());
	rpropertylist(8,state.gethv());		rpropertylist(9,state.gethl());		rpropertylist(10,state.getCpl());
	rpropertylist(11,state.getrhol());	rpropertylist(12,state.getrhov());	rpropertylist(13,state.getkl());
	rpropertylist(14,state.getmul());	rpropertylist(15,vel);				//rpropertylist(16,Qflux);
	rpropertylist(17,state.getx());		rpropertylist(18,state.getmuv());	rpropertylist(19,state.getCpv());
	rpropertylist(20,state.getkv());

	return fnRefHTCoeff(rpropertylist,Tt);

}//MBShellTubeCondenser::rhtcoeff

double	MBShellTubeCondenser::fnRefHTCoeff(Vector& rpropertylist,double Tt)
{
		double Re, Pr, Nu, h, vel, tubeod,qual, qualmin=0.05,qualmax=0.95;
		double Tsat, hfg, g=9.81,Cpl,rhol,rhov,kl,mul;
		double hfg1, den, num, Ja;
		double mur,Cpr,kr,rhor;
		double muv,Cpv,kv;
		int phase;
		double h2phase, hsuper;

		tubeod	= VGeometry(2);
		Tsat	= rpropertylist(1);		mur		= rpropertylist(3);		Cpr		= rpropertylist(4);
		kr		= rpropertylist(5);		rhor	= rpropertylist(6);		phase	= (int)rpropertylist(7);
		hfg		= rpropertylist(8)-rpropertylist(9);
		Cpl		= rpropertylist(10);	rhol	= rpropertylist(11);	rhov	= rpropertylist(12);
		kl		= rpropertylist(13);	mul		= rpropertylist(14);	vel		= rpropertylist(15);
		qual	= rpropertylist(17);	muv		= rpropertylist(18);	Cpv		= rpropertylist(19);
		kv		= rpropertylist(20);

		switch(phase)
		{
		case SUBCOOLED:
			Pr	=	(mur*Cpr*1e3)/kr;
			Re	=	(rhor*vel*tubeod)/mur;
			Nu	=	1.13*0.5*pow(Re,0.5)*Pr;
			h	=	(Nu*kr/tubeod)*1e-3;		//kW/m^2-K
			break;
		case TWOPHASE:
			double	hsubcooled;
			Ja		=	Cpl*fabs(Tsat-Tt)/hfg;
			hfg1	=	hfg*(1+0.68*Ja)*1e3;
			num		=	g*rhol*(rhol-rhov)*pow(kl,3)*hfg1;
			den		=	mul*fabs(Tsat-Tt)*tubeod;	
			if(den>0.0&&num>0.0)
				h2phase	=	7.5*0.729*pow(num/den,0.25)*1e-3;	//kW/m^2-K
			if(qual<qualmin)
			{
				double	Pr0, Re0, Nu0;
				Pr0	=	(mul*Cpl*1e3)/kl;
				Re0	=	(rhol*vel*tubeod)/mul;
				Nu0	=	1.13*0.5*pow(Re0,0.5)*Pr0;
				hsubcooled	=	(Nu0*kl/tubeod)*1e-3;
				h	=	hsubcooled + (qual/qualmin)*(h2phase-hsubcooled);
			}
			else if(qual > qualmax)
			{
				double Pr1, Re1, Nu1;
				Pr1	=	(muv*Cpv*1e3)/kr;
				Re1	=	(rhov*vel*tubeod)/muv;
				Nu1	=	1.13*pow(Re1,0.5)*Pr1;
				hsuper=	2*(Nu1*kv/tubeod)*1e-3;
				h	=	hsuper + ((1-qual)/(1-qualmax))*(h2phase-hsuper);
			}
			else
				h	=	h2phase;
			break;
		case SUPERHEATED:
			Pr	=	(mur*Cpr*1e3)/kr;
			Re	=	(rhor*vel*tubeod)/mur;
			Nu	=	1.13*pow(Re,0.5)*Pr;
			h	=	2*(Nu*kr/tubeod)*1e-3;
			break;
		}//switch phase

		h	=	3.0*h;
		return h;
}//MBShellTubeCondenser::fnRefHTCoeff

bool	MBShellTubeCondenser::fnMeanPropsTP(double P, double x1, double x2, double& rhobar, double& rhoubar)
{
	r134astate	meanstate;
	double rhol, rhov, a, b;
	double ul, uv;
	
	if(!meanstate.setstateP(P))
	{
		errorlog.Add("MBShellTubeCondenser::fnMeanPropsTP: ","Could not set state.");
		return false;
	}
	rhol = meanstate.getrhol();
	rhov = meanstate.getrhov();
	ul	 = meanstate.getul();
	uv	 = meanstate.getuv();

	b = rhov/(rhol - rhov);
	a = rhol*b;

	rhobar = (a/(x2 - x1))*log((x2 + b)/(x1 + b));
	rhoubar = a*(uv - ul) + (a/(x2 - x1))*(ul - b*(uv - ul))*log((x2 + b)/(x1 + b));

	return true;

}//MBShellTubeCondenser::fnMeanDensityTP

/*
bool	MBShellTubeCondenser::fnMeanPropsSH(double P, double h1, double h2, double dL, double& rhobar, double& rhoubar)
{
	r134astate	meanstate;
	double n = 10.0;
	double h, dh, dz, rho, rhodz = 0.0;
	double u, rhou, rhoudz = 0.0;
	dz = dL/n;
	dh = (h2 - h1)/n;
	for(int i=1;i<=n;i++)
	{
		h = h1 + (2*i-1)*dh/2.0;
		if(!meanstate.setstate(P,h))
		{
			errorlog.Add("MBShellTubeCondenser::fnMeanPropsSH: ","Could not set state.");
			return false;
		}
		rho = meanstate.getrho();
		rhodz += rho*dz;
		u	= meanstate.getu();
		rhou = rho*u;
		rhoudz += rhou;
	}

	rhobar = rhodz/dL;
	rhoubar= rhoudz/dL;

	return true;

}//MBShellTubeCondenser::fnMeanDensitySH

bool	MBShellTubeCondenser::fnMeanPropsSC(double P, double h1, double h2, double dL, double& rhobar, double& rhoubar)
{
	r134astate	meanstate;
	double n = 10.0;
	double h, dh, dz, rho, rhodz = 0.0;
	double u, rhou, rhoudz = 0.0;
	dz = dL/n;
	dh = (h2 - h1)/n;
	for(int i=1;i<=n;i++)
	{
		h = h1 + (2*i-1)*dh/2.0;
		if(!meanstate.setstate(P,h))
		{
			errorlog.Add("MBShellTubeCondenser::fnMeanPropsSC: ","Could not set state.");
			return false;
		}
		rho = meanstate.getrho();
		rhodz += rho*dz;
		u	= meanstate.getu();
		rhou = rho*u;
		rhoudz += rhou;
	}

	rhobar = rhodz/dL;
	rhoubar= rhoudz/dL;

	return true;

}//MBShellTubeCondenser::fnMeanDensitySC
*/
bool	MBShellTubeCondenser::fnMeanPropsSH(double P, double h1, double h2, double dL, double& rhobar, double& rhoubar)
{
	r134astate	meanstate;
	double rho1, rho2, u1, u2;

	if(!meanstate.setstate(P,h1))
	{
		errorlog.Add("MBShellTubeCondenser::fnMeanPropsSH: ","Could not set state.");
		return false;
	}
	rho1=	meanstate.getrho();
	u1	=	meanstate.getu();
	if(!meanstate.setstate(P,h2))
	{
		errorlog.Add("MBShellTubeCondenser::fnMeanPropsSH: ","Could not set state.");
		return false;
	}
	rho2=	meanstate.getrho();
	u2	=	meanstate.getu();

	rhobar = (rho1 + rho2)/2.0;
	rhoubar= (rho1*u1 + rho2*u2)/2.0;

	return true;

}//MBShellTubeCondenser::fnMeanDensitySH

bool	MBShellTubeCondenser::fnMeanPropsSC(double P, double h1, double h2, double dL, double& rhobar, double& rhoubar)
{
	r134astate	meanstate;
	double rho1, rho2, u1, u2;

	if(!meanstate.setstate(P,h1))
	{
		errorlog.Add("MBShellTubeCondenser::fnMeanPropsSH: ","Could not set state.");
		return false;
	}
	rho1=	meanstate.getrho();
	u1	=	meanstate.getu();
	if(!meanstate.setstate(P,h2))
	{
		errorlog.Add("MBShellTubeCondenser::fnMeanPropsSH: ","Could not set state.");
		return false;
	}
	rho2=	meanstate.getrho();
	u2	=	meanstate.getu();

	rhobar = (rho1 + rho2)/2.0;
	rhoubar= (rho1*u1 + rho2*u2)/2.0;

	return true;

}//MBShellTubeCondenser::fnMeanDensitySC

double	MBShellTubeCondenser::fnDhindt(double hrin)
{
	double dhindt;
	int step;
	step = int((dTIME*10.0))/10;
	int depth = HINDEPTH;
	if(step < depth)
	{
		dHRINOLD(step+1,hrin);
		dhindt = 0.0;
	}
	else
	{
		for(int i=1;i<depth;i++)
			dHRINOLD(i,dHRINOLD(i+1));
		dHRINOLD(depth,hrin);
		dhindt = (hrin - dHRINOLD(1))/depth;
	}

	return dhindt;

}//MBShellTubeCondenser::fnDhindt

void	MBShellTubeCondenser::fnComputeRefMassEngy(double& mass, double& engy)
{
	double hv, hl, hout, rhol;
	double rho1, rho2, rho3;
	double rhou1, rhou2, rhou3;
	double h1, h2, P, dL, L1, L2, L;
	double x1, x2, farea;
	r134astate	nstate;
	P = VState(1);
	hout = VState(4);
	L1= VState(2);
	L2= VState(3);
	L = VGeometry(3);
	farea = VGeometry(4);
	nstate.setstateP(P);
	hv = nstate.gethv();
	hl = nstate.gethl();
	rhol = nstate.getrhol();
	r134astate	meanstate;

	switch(iOPMODE)
	{
	case	SH:
		if(!meanstate.setstate(P,hout))
		{
			errorlog.Add("MBShellTubeCondenser::fnComputeRefMassEngy: ","Could not set state.");
			return;
		}
		rho1	= meanstate.getrho();
		rhou1	= rho1*meanstate.getu();
//		h1	=	VINPUTS(4);
//		h2	=	hout;
		dL	=	L;
//		fnMeanPropsSH(P,h1,h2,dL,rho1,rhou1);
		mass = rho1*dL*farea;
		engy = rhou1*dL*farea;
		break;
	case	SHTP:
		h1	=	VINPUTS(4);
		h2	=	hv;
		dL	=	L1;
		fnMeanPropsSH(P,h1,h2,dL,rho1,rhou1);
		x1	=	1.0;
		x2	=	(hout - hl)/(hv - hl);
		fnMeanPropsTP(P,x1,x2,rho2,rhou2);
		mass = farea*(rho1*L1 + rho2*(L - L1));
		engy = farea*(rhou1*L1 + rhou2*(L - L1));
		break;
	case	SHTPSC:
		h1	=	VINPUTS(4);
		h2	=	hv;
		dL	=	L1;
		fnMeanPropsSH(P,h1,h2,dL,rho1,rhou1);
		x1	=	1.0;
		x2	=	0.0;
		fnMeanPropsTP(P,x1,x2,rho2,rhou2);
		h1	=	hl;
		h2	=	hout;
		fnMeanPropsSC(P,h1,h2,L - L2,rho3,rhou3);
		mass = farea*(rho1*L1 + rho2*(L2 - L1) + rho3*(L - L2));
		engy = farea*(rhou1*L1 + rhou2*(L2 - L1) + rhou3*(L - L2));
		break;
	}

	return;
}//MBShellTubeCondenser::fnComputeRefMass

bool	MBShellTubeCondenser::fnSHDerivatives(Vector& BCs,Vector& State,Vector& Derivs,double& Qr, double& Qw)
{
	double	tubeid, tubeod, L, farea, Ahtout, Ahtin, TubeCap, NTUBES, WatCap;
	tubeid	= VGeometry(1);		tubeod	= VGeometry(2);		L		= VGeometry(3);
	farea	= VGeometry(4);		Ahtout	= VGeometry(5);		Ahtin	= VGeometry(6);
	TubeCap	= VGeometry(7);		NTUBES	= VGeometry(8);		WatCap	= VGeometry(9);

	//Working variables...heat transfer
	double	alphar, alphaw;
	//Coefficients
	double a11, a12, a21, a22, a33, a44;
	double c1, c2, c3, c4;
	//Working variables...properties
	double	rhobar, hbar, Tbar;
	double	drhodP, drhodh;
	r134astate	zone;
	//Working variables...others
	double velr;
	double velw, AwaterX, rhow, Cpw;
	rhow	= 995.0;	//kg/m^3
	Cpw		= 4.1868;	//kJ/kg-K
	AwaterX = NTUBES*MYPI*tubeid*tubeid/4.0;

	double	mwat, Twin, mrin, hrin, mrout, dhindt;
	mwat = BCs(1);	Twin = BCs(2);	mrin = BCs(3);		hrin = BCs(4);	mrout = BCs(5);	dhindt = BCs(6);
	velw	=	mwat/(rhow*AwaterX);

	double	P, hout, Tt, Tw;
	P = State(1);		hout = State(4);
	Tt	= State(7);		Tw = State(10);

	double dPdt, dhoutdt, dTtdt, dTwdt;

	hbar = (hrin + hout)/2.0;
	if(!zone.setstate(P,hbar))
	{
		errorlog.Add("MBShellTubeCondenser::fnSHDerivatives: ","Could not set state.");
		return false;
	}
	rhobar = zone.getrho();
	Tbar = zone.getT();
	drhodh = zone.getdrdh();
	drhodP = zone.getdrdP();
	velr = (fabs(mrin + mrout)/2.0)/(rhobar*farea);

	alphar = rhtcoeff(P,hbar,velr,Tt);
	alphaw = whtcoeff(Tw,velw);
	Qr = alphar*Ahtout*L*(Tbar - Tt);
	Qw = alphaw*Ahtin*L*(Tt - Tw);

	a11 = farea*L*drhodP;
	a12 = 0.5*farea*L*drhodh;
	a21 = farea*L*(hbar*drhodP - 1.0);
	a22 = 0.5*farea*L*(rhobar + hbar*drhodh);
	a33 = TubeCap*L;
	a44 = WatCap*L;

	c1 = mrin - mrout - 0.5*farea*L*drhodh*dhindt;
	c2 = mrin*hrin - mrout*hout - Qr - 0.5*farea*L*(rhobar + hbar*drhodh)*dhindt;
	c3 = Qr - Qw;
	c4 = Qw + mwat*Cpw*(Twin - Tw);

	dPdt = ((c1/a12) - (c2/a22))/((a11/a12)-(a21/a22));
	dhoutdt = ((c1/a11) - (c2/a21))/((a12/a11) - (a22/a21));
	dTtdt = c3/a33;
	dTwdt = c4/a44;

	Derivs(1,dPdt);
	Derivs(4,dhoutdt);
	Derivs(7,dTtdt);
	Derivs(10,dTwdt);

	return true;
}//MBShellTubeCondenser::fnSHDerivatives
/*
bool	MBShellTubeCondenser::fnSHDerivatives(Vector& BCs,Vector& State,Vector& Derivs,double& Qr, double& Qw)
{
	double	tubeid, tubeod, L, farea, Ahtout, Ahtin, TubeCap, NTUBES, WatCap;
	tubeid	= VGeometry(1);		tubeod	= VGeometry(2);		L		= VGeometry(3);
	farea	= VGeometry(4);		Ahtout	= VGeometry(5);		Ahtin	= VGeometry(6);
	TubeCap	= VGeometry(7);		NTUBES	= VGeometry(8);		WatCap	= VGeometry(9);

	//Working variables...heat transfer
	double	alphar, alphaw;
	//Coefficients
	double a11, a12, a21, a22, a33, a44;
	double c1, c2, c3, c4;
	//Working variables...properties
	double	rhobar, hbar, Tbar;
	double	drhodP, drhodh;
	r134astate	zone;
	//Working variables...others
	double velr;
	double velw, AwaterX, rhow, Cpw;
	rhow	= 995.0;	//kg/m^3
	Cpw		= 4.1868;	//kJ/kg-K
	AwaterX = NTUBES*MYPI*tubeid*tubeid/4.0;

	double	mwat, Twin, mrin, hrin, mrout, dhindt;
	mwat = BCs(1);	Twin = BCs(2);	mrin = BCs(3);		hrin = BCs(4);	mrout = BCs(5);	dhindt = BCs(6);
	velw	=	mwat/(rhow*AwaterX);

	double	P, hout, Tt, Tw;
	P = State(1);		hout = State(4);
	Tt	= State(7);		Tw = State(10);

	double dPdt, dhoutdt, dTtdt, dTwdt;

	hbar = hout;
	if(!zone.setstate(P,hbar))
	{
		errorlog.Add("MBShellTubeCondenser::fnSHDerivatives: ","Could not set state.");
		return false;
	}
	rhobar = zone.getrho();
	Tbar = zone.getT();
	drhodh = zone.getdrdh();
	drhodP = zone.getdrdP();
	velr = (fabs(mrin + mrout)/2.0)/(rhobar*farea);

	alphar = rhtcoeff(P,hbar,velr,Tt);
	alphaw = whtcoeff(Tw,velw);
	Qr = alphar*Ahtout*L*(Tbar - Tt);
	Qw = alphaw*Ahtin*L*(Tt - Tw);

	a11 = farea*L*drhodP;
	a12 = farea*L*drhodh;
	a21 = farea*L*(hbar*drhodP - 1.0);
	a22 = farea*L*(rhobar + hbar*drhodh);
	a33 = TubeCap*L;
	a44 = WatCap*L;

	c1 = mrin - mrout;
	c2 = mrin*hrin - mrout*hout - Qr;
	c3 = Qr - Qw;
	c4 = Qw + mwat*Cpw*(Twin - Tw);

	dPdt = ((c1/a12) - (c2/a22))/((a11/a12)-(a21/a22));
	dhoutdt = ((c1/a11) - (c2/a21))/((a12/a11) - (a22/a21));
	dTtdt = c3/a33;
	dTwdt = c4/a44;

	Derivs(1,dPdt);
	Derivs(4,dhoutdt);
	Derivs(7,dTtdt);
	Derivs(10,dTwdt);

	return true;
}//MBShellTubeCondenser::fnSHDerivatives

*/

bool	MBShellTubeCondenser::fnTPDerivatives(Vector& BCs,Vector& State,Vector& Derivs,double& Qr, double& Qw)
{
	double	tubeid, tubeod, L, farea, Ahtout, Ahtin, TubeCap, NTUBES, WatCap;
	tubeid	= VGeometry(1);		tubeod	= VGeometry(2);		L		= VGeometry(3);
	farea	= VGeometry(4);		Ahtout	= VGeometry(5);		Ahtin	= VGeometry(6);
	TubeCap	= VGeometry(7);		NTUBES	= VGeometry(8);		WatCap	= VGeometry(9);

	//Working variables...heat transfer
	double	alphar, alphaw;
	//Coefficients
	double a11, a12, a21, a22, a33, a44;
	double c1, c2, c3, c4;
	//Working variables...properties
	double	rhobar, hbar, Tbar;
	double	drhodP, drhodh;
	r134astate	zone;
	//Working variables...others
	double velr;
	double velw, AwaterX, rhow, Cpw;
	rhow	= 995.0;	//kg/m^3
	Cpw		= 4.1868;	//kJ/kg-K
	AwaterX = NTUBES*MYPI*tubeid*tubeid/4.0;

	double	mwat, Twin, mrin, hrin, mrout, dhindt;
	mwat = BCs(1);	Twin = BCs(2);	mrin = BCs(3);		hrin = BCs(4);	mrout = BCs(5);	dhindt = BCs(6);
	velw	=	mwat/(rhow*AwaterX);

	double	P, hout, Tt, Tw;
	P = State(1);		hout = State(4);
	Tt	= State(7);		Tw = State(10);

	double dPdt, dhoutdt, dTtdt, dTwdt;

	hbar = (hrin + hout)/2.0;
	if(!zone.setstate(P,hbar))
	{
		errorlog.Add("MBShellTubeCondenser::fnTPDerivatives: ","Could not set state.");
		return false;
	}
	rhobar = zone.getrho();
	Tbar = zone.getT();
	drhodh = zone.getdrdh();
	drhodP = zone.getdrdP();
	velr = (fabs(mrin + mrout)/2.0)/(rhobar*farea);

	alphar = rhtcoeff(P,hbar,velr,Tt);
	alphaw = whtcoeff(Tw,velw);
	Qr = alphar*Ahtout*L*(Tbar - Tt);
	Qw = alphaw*Ahtin*L*(Tt - Tw);

	a11 = farea*L*drhodP;
	a12 = 0.5*farea*L*drhodh;
	a21 = farea*L*(hbar*drhodP - 1.0);
	a22 = 0.5*farea*L*(rhobar + hbar*drhodh);
	a33 = TubeCap*L;
	a44 = WatCap*L;

	c1 = mrin - mrout - 0.5*farea*L*drhodh*dhindt;
	c2 = mrin*hrin - mrout*hout - Qr - 0.5*farea*L*(rhobar + hbar*drhodh)*dhindt;
	c3 = Qr - Qw;
	c4 = Qw + mwat*Cpw*(Twin - Tw);

	dPdt = ((c1/a12) - (c2/a22))/((a11/a12)-(a21/a22));
	dhoutdt = ((c1/a11) - (c2/a21))/((a12/a11) - (a22/a21));
	dTtdt = c3/a33;
	dTwdt = c4/a44;

	Derivs(1,dPdt);
	Derivs(4,dhoutdt);
	Derivs(7,dTtdt);
	Derivs(10,dTwdt);

	return true;
}//MBShellTubeCondenser::fnTPDerivatives

bool	MBShellTubeCondenser::fnSHTPDerivatives(Vector& BCs,Vector& State,Vector& Derivs,double& Qr, double& Qw)
{
	double	tubeid, tubeod, L, farea, Ahtout, Ahtin, TubeCap, NTUBES, WatCap;
	tubeid	= VGeometry(1);		tubeod	= VGeometry(2);		L		= VGeometry(3);
	farea	= VGeometry(4);		Ahtout	= VGeometry(5);		Ahtin	= VGeometry(6);
	TubeCap	= VGeometry(7);		NTUBES	= VGeometry(8);		WatCap	= VGeometry(9);

	//Working variables...heat transfer
	double	alphar1, alphar2, alphaw1, alphaw2;
	double	Qr1, Qr2, Qw1, Qw2;
	//Coefficients
	double	a11, a12, a21, a22, a31, a32, a33, a41, a42, a43;
	double	b11, b12, b13, b21, b22, b23, b31, b32, b33;
	double	d1, d2, d3;
	//Working variables...properties
	double	rho1bar, rho2bar, h1bar, h2bar, T1bar, T2bar;
	double	drhodP1, drhodh1, drhodP2, drhodh2, dhvdP;
	double	hv, rhov;
	r134astate	zone1, zone2;
	//Working variables...others
	double velr1, velr2;
	double velw, AwaterX, rhow, Cpw;
	rhow	= 995.0;	//kg/m^3
	Cpw		= 4.1868;	//kJ/kg-K
	AwaterX = NTUBES*MYPI*tubeid*tubeid/4.0;
	double	gamma;
	gamma = 0.2;

	double	mwat, Twin, mrin, hrin, mrout, dhindt;
	mwat = BCs(1);	Twin = BCs(2);	mrin = BCs(3);		hrin = BCs(4);	mrout = BCs(5);	dhindt = BCs(6);
	velw	=	mwat/(rhow*AwaterX);

	double	P, L1, hout, mdotA, Tt1, Tt2, Tw1, Tw2;
	P	= State(1);		L1	=	State(2);	hout = State(4);	mdotA	=	State(5);
	Tt1	= State(7);		Tt2	=	State(8);
	Tw1 = State(10);	Tw2	=	State(11);

	double dPdt, dhoutdt, dL1dt, dTt1dt, dTt2dt, dTw1dt, dTw2dt, detA;

	if(!zone1.setstateP(P))
	{
		errorlog.Add("MBShellTubeCondenser::fnSHTPDerivatives: ","Could not set saturated properties.");
		return false;
	}

	hv		=	zone1.gethv();
	rhov	=	zone1.getrhov();

	h1bar = (hrin + hv)/2.0;
	h2bar = (hv + hout)/2.0;

	if(!zone1.setstate(P,h1bar))
	{
		errorlog.Add("MBShellTubeCondenser::fnSHTPDerivatives: ","Could not set state in zone 1.");
		return false;
	}
	if(!zone2.setstate(P,h2bar))
	{
		errorlog.Add("MBShellTubeCondenser::fnSHTPDerivatives: ","Could not set state in zone 2.");
		return false;
	}

	rho1bar	=	zone1.getrho();
	rho2bar	=	zone2.getrho();

	T1bar = zone1.getT();
	T2bar = zone2.getT();

	drhodh1	=	zone1.getdrdh();	drhodP1	=	zone1.getdrdP();
	drhodh2	=	zone2.getdrdh();	drhodP2	=	zone2.getdrdP();
	dhvdP	=	zone1.getdhvdP();

	velr1 = (fabs(mrin + mdotA)/2.0)/(rho1bar*farea);
	velr2 = (fabs(mdotA + mrout)/2.0)/(rho2bar*farea);

	alphar1 =	rhtcoeff(P,h1bar,velr1,Tt1);
	alphar2	=	rhtcoeff(P,h2bar,velr2,Tt2);

	alphaw1	=	whtcoeff(Tw1,velw);
	alphaw2	=	whtcoeff(Tw2,velw);

	Qr1 =	alphar1*Ahtout*L1*(T1bar - Tt1);
	Qr2	=	alphar2*Ahtout*(L - L1)*(T2bar - Tt2);
	Qr	=	Qr1 + Qr2;
	Qw1 =	alphaw1*Ahtin*L1*(Tt1 - Tw1);
	Qw2	=	alphaw2*Ahtin*(L - L1)*(Tt2 - Tw2);
	Qw	=	Qw1 + Qw2;

	a11	=	farea*L1*(0.5*drhodh1*dhvdP + drhodP1);
	a12	=	farea*(rho1bar - rhov);
	a21	=	farea*L1*(0.5*rho1bar*dhvdP + 0.5*h1bar*drhodh1*dhvdP + h1bar*drhodP1 - 1.0);
	a22	=	farea*(rho1bar*h1bar - rhov*hv);
	a31	=	farea*(L - L1)*(0.5*drhodh2*dhvdP + drhodP2);
	a32	=	-farea*(rho2bar - rhov);
	a33	=	0.5*farea*(L - L1)*drhodh2;
	a41	=	farea*(L - L1)*(0.5*rho2bar*dhvdP + 0.5*h2bar*drhodh2*dhvdP + h2bar*drhodP2 - 1.0);
	a42	=	-farea*(rho2bar*h2bar - rhov*hv);
	a43	=	0.5*farea*(L - L1)*(rho2bar + h2bar*drhodh2);

	b11	=	a11 + a31;		b12 =	a12 + a32;		b13 =	a33;
	b21	=	a21 + hv*a31;	b22	=	a22 + hv*a32;	b23	=	hv*a33;
	b31	=	a41 - hv*a31;	b32	=	a42 - hv*a32;	b33	=	a43 - hv*a33;

	d1	=	mrin - mrout - 0.5*farea*L1*drhodh1*dhindt;
	d2	=	mrin*hrin - mrout*hv - Qr1 - 0.5*farea*L1*(rho1bar + h1bar*drhodh1)*dhindt;
	d3	=	mrout*(hv - hout) - Qr2;

	detA=	b11*b22*b33-b11*b23*b32-b21*b12*b33+b21*b13*b32+b31*b12*b23-b31*b13*b22;

	dPdt	=	((b22*b33-b23*b32)*d1+(-b12*b33+b13*b32)*d2+(b12*b23-b13*b22)*d3)/detA;
	dL1dt	=	((-b21*b33+b23*b31)*d1+(b11*b33-b13*b31)*d2+(-b11*b23+b13*b21)*d3)/detA;
	dhoutdt	=	((b21*b32-b22*b31)*d1+(-b11*b32+b12*b31)*d2+(b11*b22-b12*b21)*d3)/detA;

	dTt1dt	=	(Qr1 - Qw1 - TubeCap*(Tt1 - Tt2))/(TubeCap*L1);
	dTt2dt	=	(Qr2 - Qw2)/(TubeCap*(L - L1));
	dTw1dt	=	(Qw1 + mwat*Cpw*(Tw2 - Tw1) - WatCap*(Tw1 - Tw2)*dL1dt)/(WatCap*L1);
	dTw2dt	=	(Qw2 + mwat*Cpw*(Twin - Tw2))/(WatCap*(L - L1));

	mdotA	=	mrout + a31*dPdt + a32*dL1dt + a33*dhoutdt;

	Derivs(1,dPdt);		Derivs(2,dL1dt);	Derivs(4,dhoutdt);
	Derivs(5,mdotA);	Derivs(7,dTt1dt);	Derivs(8,dTt2dt);
	Derivs(10,dTw1dt);	Derivs(11,dTw2dt);

	return true;
}//MBShellTubeCondenser::fnSHTPDerivatives

bool	MBShellTubeCondenser::fnTPSCDerivatives(Vector& BCs,Vector& State,Vector& Derivs,double& Qr, double& Qw)
{
	double	tubeid, tubeod, L, farea, Ahtout, Ahtin, TubeCap, NTUBES, WatCap;
	tubeid	= VGeometry(1);		tubeod	= VGeometry(2);		L		= VGeometry(3);
	farea	= VGeometry(4);		Ahtout	= VGeometry(5);		Ahtin	= VGeometry(6);
	TubeCap	= VGeometry(7);		NTUBES	= VGeometry(8);		WatCap	= VGeometry(9);

	double detA;

	//Working variables...heat transfer
	double	alphar1, alphar2, alphaw1, alphaw2;
	double	Qr1, Qr2, Qw1, Qw2;
	//Coefficients
	double	a11, a12, a21, a22, a31, a32, a33, a41, a42, a43;
	double	b11, b12, b13, b21, b22, b23, b31, b32, b33;
	double	d1, d2, d3;
	//Working variables...properties
	double	rho1bar, rho2bar, h1bar, h2bar, T1bar, T2bar;
	double	drhodP1, drhodh1, drhodP2, drhodh2, dhldP;
	double	hl, rhol;
	r134astate	zone1, zone2;
	//Working variables...others
	double velr1, velr2;
	double velw, AwaterX, rhow, Cpw;
	rhow	= 995.0;	//kg/m^3
	Cpw		= 4.1868;	//kJ/kg-K
	AwaterX = NTUBES*MYPI*tubeid*tubeid/4.0;

	double	mwat, Twin, mrin, hrin, mrout, dhindt;
	mwat = BCs(1);	Twin = BCs(2);	mrin = BCs(3);		hrin = BCs(4);	mrout = BCs(5);	dhindt = BCs(6);
	velw	=	mwat/(rhow*AwaterX);

	double	P, L1, hout, mdotA, Tt1, Tt2, Tw1, Tw2;
	P	= State(1);		L1	=	State(2);	hout = State(4);	mdotA	=	State(5);
	Tt1	= State(7);		Tt2	=	State(8);
	Tw1 = State(10);	Tw2	=	State(11);

	double dPdt, dhoutdt, dL1dt, dTt1dt, dTt2dt, dTw1dt, dTw2dt;

	if(!zone1.setstateP(P))
	{
		errorlog.Add("MBShellTubeCondenser::fnSHTPDerivatives: ","Could not set saturated properties.");
		return false;
	}

	hl		=	zone1.gethl();
	rhol	=	zone1.getrhol();

	h1bar = (hrin + hl)/2.0;
	h2bar = (hl + hout)/2.0;

	if(!zone1.setstate(P,h1bar))
	{
		errorlog.Add("MBShellTubeCondenser::fnSHTPDerivatives: ","Could not set state in zone 1.");
		return false;
	}
	if(!zone2.setstate(P,h2bar))
	{
		errorlog.Add("MBShellTubeCondenser::fnSHTPDerivatives: ","Could not set state in zone 2.");
		return false;
	}

	rho1bar	=	zone1.getrho();
	rho2bar	=	zone2.getrho();

	T1bar = zone1.getT();
	T2bar = zone2.getT();

	drhodh1	=	zone1.getdrdh();	drhodP1	=	zone1.getdrdP();
	drhodh2	=	zone2.getdrdh();	drhodP2	=	zone2.getdrdP();
	dhldP	=	zone1.getdhldP();

	velr1 = (fabs(mrin + mdotA)/2.0)/(rho1bar*farea);
	velr2 = (fabs(mdotA + mrout)/2.0)/(rho2bar*farea);

	alphar1 =	rhtcoeff(P,h1bar,velr1,Tt1);
	alphar2	=	rhtcoeff(P,h2bar,velr2,Tt2);

	alphaw1	=	whtcoeff(Tw1,velw);
	alphaw2	=	whtcoeff(Tw2,velw);

	Qr1 =	alphar1*Ahtout*L1*(T1bar - Tt1);
	Qr2	=	alphar2*Ahtout*(L - L1)*(T2bar - Tt2);
	Qr	=	Qr1 + Qr2;
	Qw1 =	alphaw1*Ahtin*L1*(Tt1 - Tw1);
	Qw2	=	alphaw2*Ahtin*(L - L1)*(Tt2 - Tw2);
	Qw	=	Qw1 + Qw2;

	a11	=	farea*L1*(0.5*drhodh1*dhldP + drhodP1);
	a12	=	farea*(rho1bar - rhol);
	a21	=	farea*L1*(0.5*rho1bar*dhldP + 0.5*h1bar*drhodh1*dhldP + h1bar*drhodP1 - 1.0);
	a22	=	farea*(rho1bar*h1bar - rhol*hl);
	a31	=	farea*(L - L1)*(0.5*drhodh2*dhldP + drhodP2);
	a32	=	-farea*(rho2bar - rhol);
	a33	=	0.5*farea*(L - L1)*drhodh2;
	a41	=	farea*(L - L1)*(0.5*rho2bar*dhldP + 0.5*h2bar*drhodh2*dhldP + h2bar*drhodP2 - 1.0);
	a42	=	-farea*(rho2bar*h2bar - rhol*hl);
	a43	=	0.5*farea*(L - L1)*(rho2bar + h2bar*drhodh2);

	b11	=	a11 + a31;		b12 =	a12 + a32;		b13 =	a33;
	b21	=	a21 + hl*a31;	b22	=	a22 + hl*a32;	b23	=	hl*a33;
	b31	=	a41 - hl*a31;	b32	=	a42 - hl*a32;	b33	=	a43 - hl*a33;

	d1	=	mrin - mrout - 0.5*farea*L1*drhodh1*dhindt;
	d2	=	mrin*hrin - mrout*hl - Qr1 - 0.5*farea*L1*(rho1bar + h1bar*drhodh1)*dhindt;
	d3	=	mrout*(hl - hout) - Qr2;

	detA=	b11*b22*b33-b11*b23*b32-b21*b12*b33+b21*b13*b32+b31*b12*b23-b31*b13*b22;

	dPdt	=	((b22*b33-b23*b32)*d1+(-b12*b33+b13*b32)*d2+(b12*b23-b13*b22)*d3)/detA;
	dL1dt	=	((-b21*b33+b23*b31)*d1+(b11*b33-b13*b31)*d2+(-b11*b23+b13*b21)*d3)/detA;
	dhoutdt	=	((b21*b32-b22*b31)*d1+(-b11*b32+b12*b31)*d2+(b11*b22-b12*b21)*d3)/detA;

	dTt1dt	=	(Qr1 - Qw1 - TubeCap*(Tt1 - Tt2))/(TubeCap*L1);
	dTt2dt	=	(Qr2 - Qw2)/(TubeCap*(L - L1));
	dTw1dt	=	(Qw1 + mwat*Cpw*(Tw2 - Tw1) - WatCap*(Tw1 - Tw2)*dL1dt)/(WatCap*L1);
	dTw2dt	=	(Qw2 + mwat*Cpw*(Twin - Tw2))/(WatCap*(L - L1));

	mdotA	=	mrout + a31*dPdt + a32*dL1dt + a33*dhoutdt;

	Derivs(1,dPdt);		Derivs(2,dL1dt);	Derivs(4,dhoutdt);
	Derivs(5,mdotA);	Derivs(7,dTt1dt);	Derivs(8,dTt2dt);
	Derivs(10,dTw1dt);	Derivs(11,dTw2dt);

	return true;
}//MBShellTubeCondenser::fnTPSCDerivatives

bool	MBShellTubeCondenser::fnSHTPSCDerivatives(Vector& BCs,Vector& State,Vector& Derivs,double& Qr, double& Qw)
{
	double	farea, Ahtout, Ahtin, L, tubeid, tubeod, TubeCap, NTUBES, WatCap;
	tubeid	= VGeometry(1);		tubeod	= VGeometry(2);		L		= VGeometry(3);
	farea	= VGeometry(4);		Ahtout	= VGeometry(5);		Ahtin	= VGeometry(6);
	TubeCap	= VGeometry(7);		NTUBES	= VGeometry(8);		WatCap	= VGeometry(9);

	//Working variables...heat transfer
	double	alphar1, alphar2, alphar3;
	double	alphaw1, alphaw2, alphaw3;
	double	Qr1, Qr2, Qr3;
	double	Qw1, Qw2, Qw3;
	//Working variables...coefficients
	double	a11, a12, a21, a22,	a31, a32, a33, a41, a42, a43,
			a51, a53, a54,	a61, a63, a64;
	double	b11, b12, b13, b14, 
			b21, b22, b23, b24,
			b31, b32, b33, b34,
			b41, b42, b43, b44;
	double	d1, d2, d3, d4;

	//Working variables...properties
	double	rho1bar, rho2bar, rho3bar;	//average densities in each zone
	double	h1bar, h2bar, h3bar;		//average enthalpies in each zone
	double	T1bar, T2bar, T3bar;
	double	rhov, hv, rhol, hl;
	double	drhodP1, drhodh1, 
			drhodP2, drhodh2,
			drhodP3, drhodh3,
			dhvdP, dhldP,
			drhovdP, drholdP;
	r134astate	zone1, zone2, zone3;
	//Working variables...others
	double	A1, A2, A3;
	double vel1, vel2, vel3;
	double velw, AwaterX, rhow, Cpw;
	rhow	= 995.0;	//kg/m^3
	Cpw		= 4.1868;	//kJ/kg-K
	AwaterX = NTUBES*MYPI*tubeid*tubeid/4.0;
	
	double	mwat, Twin, mrin, hrin, mrout, dhindt;
	mwat = BCs(1);	Twin = BCs(2);	mrin = BCs(3);		hrin = BCs(4);	mrout = BCs(5);	dhindt = BCs(6);
	velw	=	mwat/(rhow*AwaterX);

	double	P, L1, L2, hout, mdotA, mdotB, Tt1, Tt2, Tt3, Tw1, Tw2, Tw3;
	P = State(1);		L1 = State(2);		L2 = State(3);		hout = State(4);
	mdotA = State(5);	mdotB = State(6);	Tt1	= State(7);		Tt2 = State(8);
	Tt3 = State(9);		Tw1 = State(10);	Tw2 = State(11);	Tw3 = State(12);

	double dPdt, dL1dt, dL2dt, dhoutdt, dTt1dt, dTt2dt,dTt3dt, dTw1dt, dTw2dt, dTw3dt;
	double detA;

		if(!zone1.setstateP(P))
		{
			errorlog.Add("Could not set saturated conditions in zone 1.");
			return false;
		}

		hv		= zone1.gethv();		hl	 = zone1.gethl();
		h1bar	= (hrin + hv)/2.0;		h2bar	= (hv + hl)/2.0;		h3bar	= (hl + hout)/2.0;
		rhov	= zone1.getrhov();		rhol = zone1.getrhol();

		if(!zone1.setstate(P,h1bar))
		{
			errorlog.Add("MBShellTubeCondenser::fnSHTPSCDerivatives: ","Could not set state in zone 1.");
			return false;
		}
		if(!zone2.setstate(P,h2bar))
		{
			errorlog.Add("MBShellTubeCondenser::fnSHTPSCDerivatives: ","Could not set state in zone 2.");
			return false;
		}
		if(!zone3.setstate(P,h3bar))
		{
			errorlog.Add("MBShellTubeCondenser::fnSHTPSCDerivatives: ","Could not set state in zone 3.");
			return false;
		}
		rho1bar = zone1.getrho();		rho2bar = zone2.getrho();		rho3bar = zone3.getrho();

		T1bar	= zone1.getT();		T2bar	= zone1.getTs();		T3bar	= zone3.getT();
		A1		= Ahtout*L1;		A2		= Ahtout*(L2 - L1);		A3		= Ahtout*(L - L2);	//m^2

		vel1	= (fabs(mrin + mdotA)/2.0)/(rho1bar*farea);
		vel2	= (fabs(mdotA + mdotB)/2.0)/(rho2bar*farea);
		vel3	= (fabs(mdotB + mrout)/2.0)/(rho3bar*farea);

		alphar1	= rhtcoeff(P,h1bar,vel1,Tt1);	alphar2	= rhtcoeff(P,h2bar,vel2,Tt2);	alphar3	= rhtcoeff(P,h3bar,vel3,Tt3);
		alphaw1	= whtcoeff(Tw1,velw);			alphaw2 = whtcoeff(Tw2,velw);			alphaw3 = whtcoeff(Tw3,velw);

		Qr1		= alphar1*A1*(T1bar - Tt1);		Qr2		= alphar2*A2*(T2bar - Tt2);		Qr3		= alphar3*A3*(T3bar - Tt3);
		Qr		=	Qr1 + Qr2 + Qr3;
		Qw1		= alphaw1*A1*tubeid*(Tt1 - Tw1)/tubeod;
		Qw2		= alphaw2*A2*tubeid*(Tt2 - Tw2)/tubeod;
		Qw3		= alphaw3*A3*tubeid*(Tt3 - Tw3)/tubeod;
		Qw		=	Qw1 + Qw2 + Qw3;

		drhodP1 = zone1.getdrdP();		drhodh1 = zone1.getdrdh();
		drhodP2	= zone2.getdrdP();		drhodh2	= zone2.getdrdh();
		drhodP3 = zone3.getdrdP();		drhodh3 = zone3.getdrdh();
		dhvdP	= zone1.getdhvdP();		dhldP	= zone1.getdhldP();
		drhovdP = zone1.getdrvdP();		drholdP = zone1.getdrldP();

		//Set coefficients
		a11	=	farea*L1*(0.5*drhodh1*dhvdP + drhodP1);
		a12	=	farea*(rho1bar - rhov);
		a21	=	farea*L1*(0.5*dhvdP*(rho1bar + h1bar*drhodh1) + h1bar*drhodP1 - 1.0);
		a22	=	farea*(rho1bar*h1bar - rhov*hv);
		a31	=	farea*(L2 - L1)*(0.5*drhodh2*(dhvdP + dhldP) + drhodP2);
		a32	=	-farea*(rho2bar - rhov);
		a33	=	farea*(rho2bar - rhol);
		a41	=	farea*(L2 - L1)*(0.5*(dhvdP + dhldP)*(rho2bar + h2bar*drhodh2) + h2bar*drhodP2 - 1.0);
		a42	=	-farea*(rho2bar*h2bar - rhov*hv);
		a43	=	farea*(rho2bar*h2bar - rhol*hl);
		a51	=	farea*(L - L2)*(0.5*drhodh3*dhldP + drhodP3);
		a53	=	-farea*(rho3bar - rhol);
		a54	=	0.5*farea*(L - L2)*drhodh3;
		a61	=	farea*(L - L2)*(0.5*dhldP*(rho3bar + h3bar*drhodh3) + h3bar*drhodP3 - 1.0);
		a63	=	-farea*(rho3bar*h3bar - rhol*hl);
		a64	=	0.5*farea*(L - L2)*(rho3bar + h3bar*drhodh3);

		b11	=	a11 + a31 + a51;		b12	=	a12 + a32;		b13	=	a33 + a53;		b14	=	a54;
		b21	=	a21 + hv*(a31 + a51);	b22	=	a22 + hv*a32;	b23	=	hv*(a33 + a53);	b24	=	hv*a54;
		b31	=	a41 - hv*(a31 + a51) + hl*a51;	b32	=	a42 - hv*a32;	b33 = a43 - hv*(a33 + a53) + a53*hl;	b34 = -(hv - hl)*a54;
		b41	=	a61 - hl*a51;			b42	=	0.0;			b43	=	a63 - hl*a53;	b44	=	a64 - hl*a54;

		d1	=	mrin - mrout - 0.5*farea*L1*drhodh1*dhindt;
		d2	=	mrin*hrin - mrout*hv - Qr1 - 0.5*farea*L1*(rho1bar + h1bar*drhodh1)*dhindt;
		d3	=	mrout*(hv - hl) - Qr2;
		d4	=	mrout*(hl - hout) - Qr3;

		detA = b11*b22*b33*b44-b11*b22*b34*b43-b11*b32*b23*b44+b11*b32*b24*b43
				+b11*b42*b23*b34-b11*b42*b24*b33-b21*b12*b33*b44+b21*b12*b34*b43
				+b21*b32*b13*b44-b21*b32*b14*b43-b21*b42*b13*b34+b21*b42*b14*b33
				+b31*b12*b23*b44-b31*b12*b24*b43-b31*b22*b13*b44+b31*b22*b14*b43
				+b31*b42*b13*b24-b31*b42*b14*b23-b41*b12*b23*b34+b41*b12*b24*b33
				+b41*b22*b13*b34-b41*b22*b14*b33-b41*b32*b13*b24+b41*b32*b14*b23;
		//Set derivatives
		dPdt	 = ((b22*b33*b44-b22*b34*b43-b32*b23*b44+b32*b24*b43+b42*b23*b34-b42*b24*b33)*d1
				+(-b12*b33*b44+b12*b34*b43+b32*b13*b44-b32*b14*b43-b42*b13*b34+b42*b14*b33)*d2
				+(b12*b23*b44-b12*b24*b43-b22*b13*b44+b22*b14*b43+b42*b13*b24-b42*b14*b23)*d3
				+(-b12*b23*b34+b12*b24*b33+b22*b13*b34-b22*b14*b33-b32*b13*b24+b32*b14*b23)*d4)/detA;

		dL1dt	=	((-b21*b33*b44+b21*b34*b43+b31*b23*b44-b31*b24*b43-b41*b23*b34+b41*b24*b33)*d1
				+(b11*b33*b44-b11*b34*b43-b31*b13*b44+b31*b14*b43+b41*b13*b34-b41*b14*b33)*d2
				+(-b11*b23*b44+b11*b24*b43+b21*b13*b44-b21*b14*b43-b41*b13*b24+b41*b14*b23)*d3
				+(b11*b23*b34-b11*b24*b33-b21*b13*b34+b21*b14*b33+b31*b13*b24-b31*b14*b23)*d4)/detA;

		dL2dt	=	((b21*b32*b44-b21*b34*b42-b31*b22*b44+b31*b24*b42+b41*b22*b34-b41*b24*b32)*d1
				+(-b11*b32*b44+b11*b34*b42+b31*b12*b44-b31*b14*b42-b41*b12*b34+b41*b14*b32)*d2
				+(b11*b22*b44-b11*b24*b42-b21*b12*b44+b21*b14*b42+b41*b12*b24-b41*b14*b22)*d3
				+(-b11*b22*b34+b11*b24*b32+b21*b12*b34-b21*b14*b32-b31*b12*b24+b31*b14*b22)*d4)/detA;

		dhoutdt	=((-b21*b32*b43+b21*b33*b42+b31*b22*b43-b31*b23*b42-b41*b22*b33+b41*b23*b32)*d1
				+(b11*b32*b43-b11*b33*b42-b31*b12*b43+b31*b13*b42+b41*b12*b33-b41*b13*b32)*d2
				+(-b11*b22*b43+b11*b23*b42+b21*b12*b43-b21*b13*b42-b41*b12*b23+b41*b13*b22)*d3
				+(b11*b22*b33-b11*b23*b32-b21*b12*b33+b21*b13*b32+b31*b12*b23-b31*b13*b22)*d4)/detA;

		mdotB	=	mrout + a51*dPdt + a53*dL2dt + a54*dhoutdt;
		mdotA	=	mdotB + a31*dPdt + a32*dL1dt + a33*dL2dt;

		dTt1dt	=	(Qr1 - Qw1 - TubeCap*(Tt1 - Tt2)*dL1dt)/(TubeCap*L1);
		dTt2dt	=	(Qr2 - Qw2)/(TubeCap*(L2 - L1));
		dTt3dt	=	(Qr3 - Qw3 + TubeCap*(Tt3 - Tt2)*dL2dt)/(TubeCap*(L - L2));

		dTw1dt	=	(Qw1 + mwat*Cpw*(Tw2 - Tw1) - WatCap*(Tw1 - Tw2)*dL1dt)/(WatCap*L1);
		dTw2dt	=	(Qw2 + mwat*Cpw*(Tw3 - Tw2))/(WatCap*(L2 - L1));
		dTw3dt	=	(Qw3 + mwat*Cpw*(Twin - Tw3) + WatCap*(Tw3 - Tw2)*dL2dt)/(WatCap*(L - L2));

		Derivs(1,dPdt);		Derivs(2,dL1dt);		Derivs(3,dL2dt);		Derivs(4,dhoutdt);
		Derivs(5,mdotA);	Derivs(6,mdotB);		Derivs(7,dTt1dt);		Derivs(8,dTt2dt);
		Derivs(9,dTt3dt);	Derivs(10,dTw1dt);		Derivs(11,dTw2dt);		Derivs(12,dTw3dt);

	return true;
}//MBShellTubeCondenser::fnSHTPSCDerivatives

bool	MBShellTubeCondenser::fnDerivatives(Vector& bcs, Vector& state, Vector& derivs,double& Qr, double& Qw)
{

	switch(iOPMODE)
	{
	case	SH:
		if(!fnSHDerivatives(bcs,state,derivs,Qr,Qw))
		{
			errorlog.Add("MBShellTubeCondenser::fnSH: ","Could not compute derivatives.");
			return false;
		}
		break;
	case	SHTP:
		if(!fnSHTPDerivatives(bcs,state,derivs,Qr,Qw))
		{
			errorlog.Add("MBShellTubeCondenser::fnSHTP: ","Could not compute derivatives.");
			return false;
		}
		break;
	case	SHTPSC:
		if(!fnSHTPSCDerivatives(bcs,state,derivs,Qr,Qw))
		{
			errorlog.Add("MBShellTubeCondenser::fnSHTPSC: ","Could not compute derivatives.");
			return false;
		}
		break;
	}//switch iOPMODE

	return true;
}//MBShellTubeCondenser::fnDerivatives

bool	MBShellTubeCondenser::fnIntegrateEXPE(Vector& bcs, Vector& state, double dt)
{
	Vector	derivs;
	derivs.make(12);
	double Qr, Qw;

	//State variables
	double	P, L1, L2, hout, mdotA, mdotB, Tt1, Tt2, Tt3, Tw1, Tw2, Tw3;
	P		= state(1);		L1		= state(2); 	L2		= state(3); 
	hout	= state(4);		mdotA	= state(5);		mdotB	= state(6);
	Tt1		= state(7);		Tt2		= state(8);		Tt3		= state(9);
	Tw1		= state(10);	Tw2		= state(11);	Tw3		= state(12);

	double dPdt, dL1dt, dL2dt, dhoutdt, dTt1dt, dTt2dt, dTt3dt, dTw1dt, dTw2dt, dTw3dt;

	//Compute derivatives
	if(!fnDerivatives(bcs,state,derivs,Qr,Qw))
	{
		errorlog.Add("MBShellTubeCondenser::fnSHTP: ","Could not compute derivatives.");
		return false;
	}
	dPdt	=	derivs(1);		dL1dt	=	derivs(2);		dL2dt	=	derivs(3);
	dhoutdt	=	derivs(4);		mdotA	=	derivs(5);		mdotB	=	derivs(6);
	dTt1dt	=	derivs(7);		dTt2dt	=	derivs(8);		dTt3dt	=	derivs(9);
	dTw1dt	=	derivs(10);		dTw2dt	=	derivs(11);		dTw3dt	=	derivs(12);

	P	+=	dPdt*dt;		L1	+=	dL1dt*dt;	L2	+=	dL2dt*dt;
	hout	+=	dhoutdt*dt;
	Tt1 +=	dTt1dt*dt;		Tt2	+=	dTt2dt*dt;	Tt3	+=	dTt3dt*dt;
	Tw1 +=	dTw1dt*dt;		Tw2	+=	dTw2dt*dt;	Tw3	+=	dTw3dt*dt;

	state(1,P);		state(2,L1);	state(3,L2);
	state(4,hout);	state(5,mdotA);	state(6,mdotB);
	state(7,Tt1);	state(8,Tt2);	state(9,Tt3);
	state(10,Tw1);	state(11,Tw2);	state(12,Tw3);

	double mrin, hrin, mrout, rmass, rengy;
	mrin = bcs(3);	mrout = bcs(5);
	hrin = bcs(4);
	fnComputeRefMassEngy(rmass,rengy);
	dDYNREFQTY	=	rmass;
	dTOTMASSIN	+=	mrin*dt;
	dTOTMASSOUT	+=	mrout*dt;
	dREFMASSIMB	=	(dDYNREFQTY - dORIGREFQTY)-(dTOTMASSIN-dTOTMASSOUT);

	dDYNREFENGY	=	rengy;
	dTOTENGYIN	+=	mrin*hrin*dt;
	dTOTENGYOUT	+=	mrout*hout*dt + Qr;
	dENGYIMB	=	(dDYNREFENGY - dORIGREFENGY) - (dTOTENGYIN - dTOTENGYOUT);

	return true;
}//MBShellTubeCondenser::fnIntegrateEXPE

bool	MBShellTubeCondenser::fnIntegrateEPC(Vector& bcs, Vector& state, double dt)
{
	Vector	derivs, tstate;
	derivs.make(12);
	tstate.make(12);
	double Qr, Qw;

	//State variables
	double	P, L1, L2, hout, mdotA, mdotB, Tt1, Tt2, Tt3, Tw1, Tw2, Tw3;
	P		= state(1);		L1		= state(2); 	L2		= state(3); 
	hout	= state(4);		mdotA	= state(5);		mdotB	= state(6);
	Tt1		= state(7);		Tt2		= state(8);		Tt3		= state(9);
	Tw1		= state(10);	Tw2		= state(11);	Tw3		= state(12);

	double dPdt, dL1dt, dL2dt, dhoutdt, dTt1dt, dTt2dt, dTt3dt, dTw1dt, dTw2dt, dTw3dt;
	double dP1, dL11, dL21, dho1, dTt11, dTt21, dTt31, dTw11, dTw21, dTw31;
	double dP2, dL12, dL22, dho2, dTt12, dTt22, dTt32, dTw12, dTw22, dTw32;

	//Compute derivatives
	if(!fnDerivatives(bcs,state,derivs,Qr,Qw))
	{
		errorlog.Add("MBShellTubeCondenser::fnSHTP: ","Could not compute derivatives.");
		return false;
	}
	dPdt	=	derivs(1);		dL1dt	=	derivs(2);		dL2dt	=	derivs(3);
	dhoutdt	=	derivs(4);		mdotA	=	derivs(5);		mdotB	=	derivs(6);
	dTt1dt	=	derivs(7);		dTt2dt	=	derivs(8);		dTt3dt	=	derivs(9);
	dTw1dt	=	derivs(10);		dTw2dt	=	derivs(11);		dTw3dt	=	derivs(12);

	dP1		=	dPdt*dt;		dL11	=	dL1dt*dt;		dL21	=	dL2dt*dt;
	dho1	=	dhoutdt*dt;
	dTt11	=	dTt1dt*dt;		dTt21	=	dTt2dt*dt;		dTt31	=	dTt3dt*dt;
	dTw11	=	dTw1dt*dt;		dTw21	=	dTw2dt*dt;		dTw31	=	dTw3dt*dt;

	tstate(1,P+dP1);		tstate(2,L1+dL11);		tstate(3,L1+dL21);
	tstate(4,hout+dho1);	tstate(5,mdotA);		tstate(6,mdotB);
	tstate(7,Tt1+dTt11);	tstate(8,Tt2+dTt21);	tstate(9,Tt3+dTt31);
	tstate(10,Tw1+dTw11);	tstate(11,Tw2+dTw21);	tstate(12,Tw3+dTw31);

	//Compute derivatives
	if(!fnDerivatives(bcs,state,derivs,Qr,Qw))
	{
		errorlog.Add("MBShellTubeCondenser::fnSHTP: ","Could not compute derivatives.");
		return false;
	}
	dPdt	=	derivs(1);		dL1dt	=	derivs(2);		dL2dt	=	derivs(3);
	dhoutdt	=	derivs(4);		mdotA	=	derivs(5);		mdotB	=	derivs(6);
	dTt1dt	=	derivs(7);		dTt2dt	=	derivs(8);		dTt3dt	=	derivs(9);
	dTw1dt	=	derivs(10);		dTw2dt	=	derivs(11);		dTw3dt	=	derivs(12);

	dP2		=	dPdt*dt;			dL12	=	dL1dt*dt;		dL22	=	dL2dt*dt;
	dho2	=	dhoutdt*dt;
	dTt12	=	dTt1dt*dt;			dTt22	=	dTt2dt*dt;		dTt32	=	dTt3dt*dt;
	dTw12	=	dTw1dt*dt;			dTw22	=	dTw2dt*dt;		dTw32	=	dTw3dt*dt;

	P		+=	(dP1+dP2)/2.0;			L1	+=	(dL11 + dL12)/2.0;		L2	+=	(dL21 + dL22)/2.0;
	hout	+=	(dho1 + dho2)/2.0;
	Tt1		+=	(dTt11 + dTt12)/2.0;	Tt2	+=	(dTt21 + dTt22)/2.0;	Tt3	+=	(dTt31 + dTt32)/2.0;
	Tw1		+=	(dTw11 + dTw12)/2.0;	Tw2	+=	(dTw21 + dTw22)/2.0;	Tw3	+=	(dTw31 + dTw32)/2.0;

	state(1,P);		state(2,L1);	state(3,L2);
	state(4,hout);	state(5,mdotA);	state(6,mdotB);
	state(7,Tt1);	state(8,Tt2);	state(9,Tt3);
	state(10,Tw1);	state(11,Tw2);	state(12,Tw3);

	double mrin, hrin, mrout, rmass, rengy;
	mrin = bcs(3);	mrout = bcs(5);
	hrin = bcs(4);
	fnComputeRefMassEngy(rmass,rengy);
	dDYNREFQTY	=	rmass;
	dTOTMASSIN	+=	mrin*dt;
	dTOTMASSOUT	+=	mrout*dt;
	dREFMASSIMB	=	(dDYNREFQTY - dORIGREFQTY)-(dTOTMASSIN-dTOTMASSOUT);

	dDYNREFENGY	=	rengy;
	dTOTENGYIN	+=	mrin*hrin*dt;
	dTOTENGYOUT	+=	mrout*hout*dt + Qr;
	dENGYIMB	=	(dDYNREFENGY - dORIGREFENGY) - (dTOTENGYIN - dTOTENGYOUT);

	return true;
}//MBShellTubeCondenser::fnIntegrateEPC

bool	MBShellTubeCondenser::fnChkModeChange(Vector& NewState, int& NewMode)
{
	r134astate	exitstate;
	double P, hout, hv, hl;
	P	=	NewState(1);
	hout=	NewState(4);
	if(!exitstate.setstateP(P))
	{
		errorlog.Add("MBShellTubeCondenser::fnChkModeChange: ","Could not set exit state.");
		return false;
	}
	hv = exitstate.gethv();		hl = exitstate.gethl();
	switch(iOPMODE)
	{
	case	SH:
		if(hout < hv)				//output is TP
		{			NewMode	=	SHTP;		iMODECHANGE = true;	iMODECHANGETYPE	=	SH_SHTP;		}
		else
		{			NewMode	=	iOPMODE;	iMODECHANGE = false;									}
		break;
	case	SHTP:
		if(hout > hv)				//output is SH
		{			NewMode	=	SH;			iMODECHANGE = true;	iMODECHANGETYPE	=	SHTP_SH;		}
		else if(hout < hl)			//output is SC
		{			NewMode	=	SHTPSC;		iMODECHANGE = true;	iMODECHANGETYPE	=	SHTP_SHTPSC;	}
		else
		{			NewMode	=	iOPMODE;	iMODECHANGE = false;									}
		break;
	case	SHTPSC:
		if(hout > hl)				//output is TP
		{			NewMode	=	SHTP;		iMODECHANGE = true;	iMODECHANGETYPE =	SHTPSC_SHTP;	}
		else
		{			NewMode	=	iOPMODE;	iMODECHANGE = false;									}
		break;
	}

	return true;
}//MBShellTubeCondenser::fnChkModeChange

bool	MBShellTubeCondenser::fnComputeDt(double P, double hout, double dPdt, double dhoutdt, double& dt)
{
	char errstr[512];

	r134astate nstate;
	double Pnew, honew, hbound, res1, res2, temp;
	bool done = false;
	double dt1 = 0.001, dt2 = 0.005;

	Pnew	= P + dPdt*dt1;	
	honew	= hout + dhoutdt*dt1;
	if(!nstate.setstate(Pnew,honew))
	{
		sprintf(errstr,"Could not set state at %f kPa, %f kJ/kg.",Pnew,honew);
		errorlog.Add("MBShelTubeCondenser::fnComputeDt: ",errstr);
		return false;
	}
	if(iMODECHANGETYPE == SH_SHTP || iMODECHANGETYPE == SHTP_SH)
		hbound = nstate.gethv();
	else
		hbound = nstate.gethl();
	res1 = fabs(honew - hbound);

	while(!done)
	{
		Pnew	= P + dPdt*dt2;
		honew	= hout + dhoutdt*dt2;
		if(!nstate.setstate(Pnew,honew))
		{
			sprintf(errstr,"Could not set state at %f kPa, %f kJ/kg.",Pnew,honew);
			errorlog.Add("MBShelTubeCondenser::fnComputeDt: ","Could not set state.");
			return false;
		}
		if(iMODECHANGETYPE == SH_SHTP || iMODECHANGETYPE == SHTP_SH)
			hbound = nstate.gethv();
		else
			hbound = nstate.gethl();
		res2 = fabs(honew - hbound);

		if(fabs(res1 - res2) < 1e-9)
			done = true;
		else
		{
			temp = dt1 - res1*(dt1 - dt2)/(res1 - res2);
			dt1	 = dt2;
			dt2	 = temp;
			res1 = res2;
		}
	}//while !done

	dt = (dt1>dt2)?dt1:dt2;

	return true;
}//MBShellTubeCondenser::fnComputeDt

bool	MBShellTubeCondenser::fnAdvance1step(double& dt, double dhindt)
{
	//Boundary conditions
	double	mwat, Twin, mrin, hrin, mrout;
	mwat = VINPUTS(1);	Twin = VINPUTS(2); mrin = VINPUTS(3); hrin = VINPUTS(4);
	mrout= VINPUTS(5);
	//State variables
	double	P, L1, L2, hout, mdotA, mdotB, Tt1, Tt2, Tt3, Tw1, Tw2, Tw3;
	P		= VState(1);	L1		= VState(2); 	L2		= VState(3); 
	hout	= VState(4);	mdotA	= VState(5);	mdotB	= VState(6);
	Tt1		= VState(7);	Tt2		= VState(8);	Tt3		= VState(9);
	Tw1		= VState(10);	Tw2		= VState(11);	Tw3		= VState(12);
	double	L;
	L = VGeometry(3);

	Vector bcs, state, derivs;
	bcs.make(6);	state.make(12);	derivs.make(12);
	bcs(1,mwat);	bcs(2,Twin);	bcs(3,mrin);	bcs(4,hrin);	bcs(5,mrout);	bcs(6,dhindt);
	state = VState;
	int	NewMode;
	double	epsl;
	epsl = (dTIME < 5.0)?EPSL0:EPSL1;

	if(!fnIntegrateEXPE(bcs,state,dt))
	{
		errorlog.Add("MBShellTubeCondenser::fnAdvance1step: ","Could not integrate.");
		return false;
	}
	if(!fnChkModeChange(state,NewMode))
		//SPEED IMPROVEMENT IDEA...pass hv, hl from Integrate to ChkModeChange
		//to avoid re-evaluation of state properties
	{
		errorlog.Add("MBShellTubeCondenser::fnAdvance1step: ","Could not check mode.");
		return false;
	}
	if(iMODECHANGE)
	{
		VState = state;
		iOPMODE	=	NewMode;
		if(iMODECHANGETYPE == SH_SHTP)
		{
			VState(2,L - epsl);	VState(5,mrout);
			VState(8,Tt1);	VState(11,Tw1);
		}
		if(iMODECHANGETYPE == SHTP_SHTPSC)
		{
			VState(3,L - epsl);	VState(6,mrout);
			VState(9,Tt2);	VState(12,Tw2);
		}
		iMODECHANGE = false;
	}//if iMODECHANGE
	else
		VState = state;
	VOUTPUTS = VState;

	return true;
}//MBShellTubeCondenser::fnAdvance1step

bool	MBShellTubeCondenser::fnAdvance1sec()
{
	char errstr[512];

	double tstart, tend, dt;
	tstart = 0.0;	tend = 1.0;
	dt = dTSTEP;

	r134astate	exitstate, instate;
	double	exitx, inx;

	double dhindt;
	dhindt = fnDhindt(VINPUTS(4));
	printf("dhindt = %f%c",dhindt,13);

	while(tstart < tend)
	{
		if(!fnAdvance1step(dt,dhindt))
		{
			sprintf(errstr,"Could not advance condenser state at time %f",tstart);
			errorlog.Add("MBShellTubeCondenser::fnAdvance1sec: ",errstr);
			return false;
		}
		instate.setstate(VState(1),VINPUTS(4));
		inx = instate.getx();
		exitstate.setstate(VState(1),VState(4));
		exitx = exitstate.getx();
		tstart	+=	dt;
		dTIME	+=	dt;
		dt = (tstart + dt)>tend?(tend - tstart):dTSTEP;
	}

	return true;
}//MBShellTubeCondenser::fnAdvance1sec

/************************************************************************/
/**********************E V A P O R A T O R*******************************/
/************************************************************************/

MBShellTubeEvaporator::MBShellTubeEvaporator()
{
	VINPUTS.make(5);
	VOUTPUTS.make(12);
	VGeometry.make(10);
	VState.make(12);
	rpropertylist.make(20);
	wpropertylist.make(5);
}//MBShellTubeEvaporator::MBShellTubeEvaporator

bool	MBShellTubeEvaporator::fnDefineGeometry(const char* fname)
{
	ifstream	infile;
	infile.open(fname);

	int N, scrap;
	double tubeid, tubeod, tubelen, tubeCp, tuberho, foulingfactor, shellid;

	infile >> iHXTYPE >> scrap >> N;
	infile >> tubeid >> tubeod >> tubelen >> tubeCp >> tuberho >> foulingfactor >> shellid;
	infile.close();

	NTUBES = N;

	double TotShellVolume, TotHTArea, TotTubeCap, TotWatCap, farea;
	double waterrho = 995.0, waterCp = 4.1868;

	farea			=	(MYPI/4)*(shellid*shellid-NTUBES*tubeod*tubeod);
	TotShellVolume	=	TotVolume = farea;		//Volume per m
	TotHTArea		=	NTUBES*MYPI*tubeod;		//Outside ht area per m
	TotTubeCap		=	NTUBES*(MYPI/4)*(tubeod*tubeod-tubeid*tubeid)*tuberho*tubeCp;	//Tube capacitance per m
	TotWatCap		=	NTUBES*(MYPI/4)*tubeid*tubeid*waterrho*waterCp;					//Water capacitance per m

	VGeometry(1,tubeid);		VGeometry(2,tubeod);			VGeometry(3,tubelen);
	VGeometry(4,farea);			VGeometry(5,TotHTArea);			VGeometry(6,TotHTArea*tubeid/tubeod);
	VGeometry(7,TotTubeCap);	VGeometry(8,double(NTUBES));	VGeometry(9,TotWatCap);
	VGeometry(10,foulingfactor);

	return true;
}//MBShellTubeEvaporator::fnDefineGeometry

bool	MBShellTubeEvaporator::fnInitialize(Vector& init)
{
	iN_INPUTS = 5;	iN_OUTPUTS = 12;

	VState(1,init(1));		//Pressure, kPa
	VState(2,init(2));		//L1, m
	VState(3,init(3));		//L2, m
	VState(4,init(4));		//hout, m
	VState(5,init(5));		//mdotA, kg/s
	VState(6,init(6));		//mdotB, kg/s	
	VState(7,init(7));		//Tt1, C
	VState(8,init(8));		//Tt2, C
	VState(9,init(9));		//Tt3, C
	VState(10,init(10));	//Tw1, C
	VState(11,init(11));	//Tw2, C
	VState(12,init(12));	//Tw3, C

	double	L1, L2, L;
	L1	=	init(2);
	L2	=	init(3);
	L = VGeometry(3);

	if(L1 < 0.0 || L2 < 0.0 || L1 > L || L2 > L)
	{
		errorlog.Add("MBShellTubeEvaporator::fnInitialize: ","Invalid initialization.");
		return false;
	}

	r134astate	initstate;
	double P, hout, hv, hl;
	P = init(1);
	hout = init(4);
	if(!initstate.setstateP(P))
	{
		errorlog.Add("MBShellTubeEvaporator::fnInitialize: ","Could  not set initial state.");
		return false;
	}

	hv = initstate.gethv();	hl = initstate.gethl();

	if(hout < hv)
		iOPMODE	=	TP;
	else
		iOPMODE	=	TPSH;

	double rho, farea, u;
	farea = VGeometry(4);

	if(!initstate.setstate(P,hout))
	{
		errorlog.Add("MBShellTubeEvaporator::fnInitialize: ","Could not set initial exit state.");
		return false;
	}
	rho	=	initstate.getrho();
	u	=	initstate.getu();

	dORIGREFQTY =	rho*L*farea;
	dORIGREFENGY=	dORIGREFQTY*u;

	dTIME	=	0.0;
	dTSTEP	=	0.01;
	dTOTMASSIN	=	0.0;
	dTOTMASSOUT	=	0.0;
	dTOTENGYIN	=	0.0;
	dTOTENGYOUT	=	0.0;
	dQflux		=	0.0;
	dHRINOLD.make(HINDEPTH);

	return true;

}//MBShellTubeEvaporator::fnInitialize

double	MBShellTubeEvaporator::whtcoeff(double Tw, double vel)
{
	watstate water;
	water.setstate(Tw);

	wpropertylist(1,water.getmu());
	wpropertylist(2,water.getk());
	wpropertylist(3,water.getCp());
	wpropertylist(4,water.getrho());
	wpropertylist(5,vel);

	return fnWatHTCoeff(wpropertylist);
}//MBShellTubeEvaporator::whtcoeff

double	MBShellTubeEvaporator::fnWatHTCoeff(Vector& wpropertylist)
{

	double	Pr,Re,Nu,rho,mu,Cp,k,h,f,vel,tubeid;
	mu	=	wpropertylist(1);
	k	=	wpropertylist(2);
	Cp	=	wpropertylist(3);
	rho	=	wpropertylist(4);
	vel	=	wpropertylist(5);

	tubeid = VGeometry(1);

	Pr	= mu*Cp*1e3/k;
	Re	= rho*vel*tubeid/mu;

	f	=	1/pow(0.79*log(Re)-1.64,2);
	Nu	=	(f*Re*Pr/8)/(1.07+12.7*sqrt(f/8)*(pow(Pr,0.67)-1));
	h	=	2.4*k*1e-3*Nu/tubeid;

	h	=	3.0*h;
	return h;
}//MBShellTubeEvaporator::fnWatHTCoeff

double	MBShellTubeEvaporator::rhtcoeff(double P, double h, double vel, double Tt)
{
	r134astate state;
	state.setstate(P,h);

	rpropertylist(1,state.getTs());		rpropertylist(3,state.getmu());		rpropertylist(4,state.getCp());
	rpropertylist(5,state.getk());		rpropertylist(6,state.getrho());	rpropertylist(7,state.getphase());
	rpropertylist(8,state.gethv());		rpropertylist(9,state.gethl());		rpropertylist(10,state.getCpl());
	rpropertylist(11,state.getrhol());	rpropertylist(12,state.getrhov());	rpropertylist(13,state.getkl());
	rpropertylist(14,state.getmul());	rpropertylist(15,vel);				rpropertylist(16,dQflux);
	rpropertylist(17,state.getx());		rpropertylist(18,state.getmuv());	rpropertylist(19,state.getCpv());
	rpropertylist(20,state.getkv());

	return fnRefHTCoeff(rpropertylist,Tt);

}//MBShellTubeEvaporator::rhtcoeff

double	MBShellTubeEvaporator::fnRefHTCoeff(Vector& rpropertylist,double Tt)
{
		double Re, Pr, Nu, h, vel, tubeod,Qflux,qual, qualmin=0.05,qualmax=0.95;
		double Tsat, hfg, g=9.81,Cpl,rhol,rhov,kl,mul;
		double mur,Cpr,kr,rhor;
		double muv,Cpv,kv;
		int phase;
		double h2phase, hsuper;

		tubeod	= VGeometry(2);
		Tsat	= rpropertylist(1);		mur		= rpropertylist(3);		Cpr		= rpropertylist(4);
		kr		= rpropertylist(5);		rhor	= rpropertylist(6);		phase	= (int)rpropertylist(7);
		hfg		= rpropertylist(8)-rpropertylist(9);
		Cpl		= rpropertylist(10);	rhol	= rpropertylist(11);	rhov	= rpropertylist(12);
		kl		= rpropertylist(13);	mul		= rpropertylist(14);	vel		= rpropertylist(15);
		Qflux	= rpropertylist(16);	qual	= rpropertylist(17);	muv		= rpropertylist(18);
		Cpv		= rpropertylist(19);	kv		= rpropertylist(20);

		switch(phase)
		{
		case SUBCOOLED:
			Pr	=	(mur*Cpr*1e3)/kr;
			Re	=	(rhor*vel*tubeod)/mur;
			Nu	=	1.13*0.5*pow(Re,0.5)*Pr;
			h	=	(Nu*kr/tubeod)*1e-3;		//kW/m^2-K
			break;
		case TWOPHASE:
				h2phase =	12.492 + Qflux/2.22;	//kW/m^2-K
				if(qual<=qualmax)
					h	=	h2phase;
				else
				{
					double Pr0, Re0, Nu0;
					Pr0	= (muv*Cpv*1e3)/kv;
					Re0	= (rhov*vel*tubeod)/muv;
					Nu0 = 1.13*0.5*pow(Re0,0.5)*Pr0;
					hsuper = 4*(Nu0*kv/tubeod)*1e-3;
					h	=	hsuper + ((1-qual)/(1-qualmax))*(h2phase-hsuper);
				}
			break;
		case SUPERHEATED:
			Pr	=	(mur*Cpr*1e3)/kr;
			Re	=	(rhor*vel*tubeod)/mur;
			Nu	=	1.13*pow(Re,0.5)*Pr;
			h	=	2*(Nu*kr/tubeod)*1e-3;
			break;
		}//switch phase

		h	=	1.25*h;
		return h;
}//MBShellTubeEvaporator::fnRefHTCoeff


bool	MBShellTubeEvaporator::fnMeanPropsTP(double P, double x1, double x2, double& rhobar, double& rhoubar)
{
	r134astate	meanstate;
	double rhol, rhov, a, b;
	double ul, uv;
	
	if(!meanstate.setstateP(P))
	{
		errorlog.Add("MBShellTubeCondenser::fnMeanPropsTP: ","Could not set state.");
		return false;
	}
	rhol = meanstate.getrhol();
	rhov = meanstate.getrhov();
	ul	 = meanstate.getul();
	uv	 = meanstate.getuv();

	b = rhov/(rhol - rhov);
	a = rhol*b;

	rhobar = (a/(x2 - x1))*log((x2 + b)/(x1 + b));
	rhoubar = a*(uv - ul) + (a/(x2 - x1))*(ul - b*(uv - ul))*log((x2 + b)/(x1 + b));

	return true;

}//MBShellTubeEvaporator::fnMeanDensityTP

bool	MBShellTubeEvaporator::fnMeanPropsSH(double P, double h1, double h2, double dL, double& rhobar, double& rhoubar)
{
	r134astate	meanstate;
	double rho1, rho2, u1, u2;

	if(!meanstate.setstate(P,h1))
	{
		errorlog.Add("MBShellTubeCondenser::fnMeanPropsSH: ","Could not set state.");
		return false;
	}
	rho1=	meanstate.getrho();
	u1	=	meanstate.getu();
	if(!meanstate.setstate(P,h2))
	{
		errorlog.Add("MBShellTubeCondenser::fnMeanPropsSH: ","Could not set state.");
		return false;
	}
	rho2=	meanstate.getrho();
	u2	=	meanstate.getu();

	rhobar = (rho1 + rho2)/2.0;
	rhoubar= (rho1*u1 + rho2*u2)/2.0;

	return true;

}//MBShellTubeEvaporator::fnMeanDensitySH

double	MBShellTubeEvaporator::fnDhindt(double hrin)
{
	double dhindt, hout;
	int step;
	step = int((dTIME*10.0))/10;
	int depth = HINDEPTH;
	hout = VState(4);
	if(step < depth)
	{
		dHRINOLD(step+1,hrin);
		dhindt = 0.0;
	}
	else
	{
		for(int i=1;i<depth;i++)
			dHRINOLD(i,dHRINOLD(i+1));
		dHRINOLD(depth,hrin);
		dhindt = (hrin - dHRINOLD(1))/depth;
	}

	return dhindt;

}//MBShellTubeEvaporator::fnDhindt
void	MBShellTubeEvaporator::fnComputeRefMassEngy(double& mass, double& engy)
{
	double hv, hl, hout, rhol;
	double rho1, rho2;
	double rhou1, rhou2;
	double h1, h2, P, dL, L1, L2, L;
	double x1, x2, farea;
	r134astate	nstate;
	P = VState(1);
	hout = VState(4);
	L1= VState(2);
	L2= VState(3);
	L = VGeometry(3);
	farea = VGeometry(4);
	nstate.setstateP(P);
	hv = nstate.gethv();
	hl = nstate.gethl();
	rhol = nstate.getrhol();
	r134astate	meanstate;

	switch(iOPMODE)
	{
	case	TP:
		if(!meanstate.setstate(P,hout))
		{
			errorlog.Add("MBShellTubeEvaporator::fnComputeRefMassEngy: ","Could not set state.");
			return;
		}
		rho1 = meanstate.getrho();
		rhou1 = rho1*meanstate.getu();
//		h1	=	VINPUTS(4);
//		h2	=	hout;
		dL	=	L;
//		x1	=	(h1 - hl)/(hv - hl);
//		x2	=	(hout - hl)/(hv - hl);
//		fnMeanPropsTP(P,x1,x2,rho1,rhou1);
		mass = rho1*dL*farea;
		engy = rhou1*dL*farea;
		break;
	case	TPSH:
		h1	=	VINPUTS(4);
		h2	=	hv;
		dL	=	L1;
		x1	=	(h1 - hl)/(hv - hl);
		x2	=	1.0;
		fnMeanPropsTP(P,x1,x2,rho1,rhou1);
		h1	=	hv;
		h2	=	hout;
		dL	=	L - L1;
		fnMeanPropsSH(P,h1,h2,dL,rho2,rhou2);
		mass = farea*(rho1*L1 + rho2*(L - L1));
		engy = farea*(rhou1*L1 + rhou2*(L - L1));
		break;
	}
	return;
}//MBShellTubeEvaporator::fnComputeRefMass
/*
bool	MBShellTubeEvaporator::fnTPDerivatives(Vector& BCs,Vector& State,Vector& Derivs,double& Qr, double& Qw)
{
	double	tubeid, tubeod, L, farea, Ahtout, Ahtin, TubeCap, NTUBES, WatCap;
	tubeid	= VGeometry(1);		tubeod	= VGeometry(2);		L		= VGeometry(3);
	farea	= VGeometry(4);		Ahtout	= VGeometry(5);		Ahtin	= VGeometry(6);
	TubeCap	= VGeometry(7);		NTUBES	= VGeometry(8);		WatCap	= VGeometry(9);

	//Working variables...heat transfer
	double	alphar, alphaw;
	//Coefficients
	double a11, a12, a21, a22, a33, a44;
	double c1, c2, c3, c4;
	//Working variables...properties
	double	rhobar, hbar, Tbar;
	double	drhodP, drhodh;
	r134astate	zone;
	//Working variables...others
	double velr;
	double velw, AwaterX, rhow, Cpw;
	rhow	= 995.0;	//kg/m^3
	Cpw		= 4.1868;	//kJ/kg-K
	AwaterX = NTUBES*MYPI*tubeid*tubeid/4.0;

	double	mwat, Twin, mrin, hrin, mrout, dhindt;
	mwat = BCs(1);	Twin = BCs(2);	mrin = BCs(3);		hrin = BCs(4);	mrout = BCs(5);	dhindt = BCs(6);
	velw	=	mwat/(rhow*AwaterX);

	double	P, hout, Tt, Tw;
	P = State(1);		hout = State(4);
	Tt	= State(7);		Tw = State(10);

	double dPdt, dhoutdt, dTtdt, dTwdt;

	hbar = (hrin + hout)/2.0;
	if(!zone.setstate(P,hbar))
	{
		errorlog.Add("MBShellTubeEvaporator::fnTPDerivatives: ","Could not set state.");
		return false;
	}
	rhobar = zone.getrho();
	Tbar = zone.getT();
	drhodh = zone.getdrdh();
	drhodP = zone.getdrdP();
	velr = (fabs(mrin + mrout)/2.0)/(rhobar*farea);

	alphar = rhtcoeff(P,hbar,velr,Tt);
	alphaw = whtcoeff(Tw,velw);
	Qr = alphar*Ahtout*L*(Tbar - Tt);
	dQflux = fabs(Qr)/Ahtout;
	Qw = alphaw*Ahtin*L*(Tt - Tw);

	a11 = farea*L*drhodP;
	a12 = 0.5*farea*L*drhodh;
	a21 = farea*L*(hbar*drhodP - 1.0);
	a22 = 0.5*farea*L*(rhobar + hbar*drhodh);
	a33 = TubeCap*L;
	a44 = WatCap*L;

	c1 = mrin - mrout - 0.5*farea*L*drhodh*dhindt;
	c2 = mrin*hrin - mrout*hout - Qr - 0.5*farea*L*(rhobar + hbar*drhodh)*dhindt;
	c3 = Qr - Qw;
	c4 = Qw + mwat*Cpw*(Twin - Tw);

	dPdt = ((c1/a12) - (c2/a22))/((a11/a12)-(a21/a22));
	dhoutdt = ((c1/a11) - (c2/a21))/((a12/a11) - (a22/a21));
	dTtdt = c3/a33;
	dTwdt = c4/a44;

	Derivs(1,dPdt);
	Derivs(4,dhoutdt);
	Derivs(7,dTtdt);
	Derivs(10,dTwdt);

	return true;
}//MBShellTubeEvaporator::fnTPDerivatives
*/
bool	MBShellTubeEvaporator::fnTPDerivatives(Vector& BCs,Vector& State,Vector& Derivs,double& Qr, double& Qw)
{
	double	tubeid, tubeod, L, farea, Ahtout, Ahtin, TubeCap, NTUBES, WatCap;
	tubeid	= VGeometry(1);		tubeod	= VGeometry(2);		L		= VGeometry(3);
	farea	= VGeometry(4);		Ahtout	= VGeometry(5);		Ahtin	= VGeometry(6);
	TubeCap	= VGeometry(7);		NTUBES	= VGeometry(8);		WatCap	= VGeometry(9);

	//Working variables...heat transfer
	double	alphar, alphaw;
	//Coefficients
	double a11, a12, a21, a22, a33, a44;
	double c1, c2, c3, c4;
	//Working variables...properties
	double	rhobar, hbar, Tbar;
	double	drhodP, drhodh;
	r134astate	zone;
	//Working variables...others
	double velr;
	double velw, AwaterX, rhow, Cpw;
	rhow	= 995.0;	//kg/m^3
	Cpw		= 4.1868;	//kJ/kg-K
	AwaterX = NTUBES*MYPI*tubeid*tubeid/4.0;

	double	mwat, Twin, mrin, hrin, mrout, dhindt;
	mwat = BCs(1);	Twin = BCs(2);	mrin = BCs(3);		hrin = BCs(4);	mrout = BCs(5);	dhindt = BCs(6);
	velw	=	mwat/(rhow*AwaterX);

	double	P, hout, Tt, Tw;
	P = State(1);		hout = State(4);
	Tt	= State(7);		Tw = State(10);

	double dPdt, dhoutdt, dTtdt, dTwdt;

	hbar = hout;
	if(!zone.setstate(P,hbar))
	{
		errorlog.Add("MBShellTubeEvaporator::fnTPDerivatives: ","Could not set state.");
		return false;
	}
	rhobar	= zone.getrho();
	Tbar	= zone.getT();
	drhodh	= zone.getdrdh();
	drhodP	= zone.getdrdP();
	velr	= (fabs(mrin + mrout)/2.0)/(rhobar*farea);

	alphar	= rhtcoeff(P,hbar,velr,Tt);
	alphaw	= whtcoeff(Tw,velw);
	Qr		= alphar*Ahtout*L*(Tbar - Tt);
	dQflux	= fabs(Qr)/Ahtout;
	Qw		= alphaw*Ahtin*L*(Tt - Tw);

	a11 = farea*L*drhodP;
	a12 = farea*L*drhodh;
	a21 = farea*L*(hbar*drhodP - 1.0);
	a22 = farea*L*(rhobar + hbar*drhodh);
	a33 = TubeCap*L;
	a44 = WatCap*L;

	c1 = mrin - mrout;
	c2 = mrin*hrin - mrout*hout - Qr;
	c3 = Qr - Qw;
	c4 = Qw + mwat*Cpw*(Twin - Tw);

	dPdt = ((c1/a12) - (c2/a22))/((a11/a12)-(a21/a22));
	dhoutdt = ((c1/a11) - (c2/a21))/((a12/a11) - (a22/a21));
	dTtdt = c3/a33;
	dTwdt = c4/a44;

	Derivs(1,dPdt);
	Derivs(4,dhoutdt);
	Derivs(7,dTtdt);
	Derivs(10,dTwdt);

	return true;
}//MBShellTubeEvaporator::fnTPDerivatives

bool	MBShellTubeEvaporator::fnTPSHDerivatives(Vector& BCs,Vector& State,Vector& Derivs,double& Qr, double& Qw)
{
	double	tubeid, tubeod, L, farea, Ahtout, Ahtin, TubeCap, NTUBES, WatCap;
	tubeid	= VGeometry(1);		tubeod	= VGeometry(2);		L		= VGeometry(3);
	farea	= VGeometry(4);		Ahtout	= VGeometry(5);		Ahtin	= VGeometry(6);
	TubeCap	= VGeometry(7);		NTUBES	= VGeometry(8);		WatCap	= VGeometry(9);

	//Working variables...heat transfer
	double	alphar1, alphar2, alphaw1, alphaw2;
	double	Qr1, Qr2, Qw1, Qw2;
	//Coefficients
	double	a11, a12, a21, a22, a31, a32, a33, a41, a42, a43;
	double	b11, b12, b13, b21, b22, b23, b31, b32, b33;
	double	d1, d2, d3;
	//Working variables...properties
	double	rho1bar, rho2bar, h1bar, h2bar, T1bar, T2bar;
	double	drhodP1, drhodh1, drhodP2, drhodh2, dhvdP;
	double	hv, rhov;
	r134astate	zone1, zone2;
	//Working variables...others
	double velr1, velr2;
	double velw, AwaterX, rhow, Cpw;
	rhow	= 995.0;	//kg/m^3
	Cpw		= 4.1868;	//kJ/kg-K
	AwaterX = NTUBES*MYPI*tubeid*tubeid/4.0;
	double	gamma;
	gamma = 0.2;

	double	mwat, Twin, mrin, hrin, mrout, dhindt;
	mwat = BCs(1);	Twin = BCs(2);	mrin = BCs(3);		hrin = BCs(4);	mrout = BCs(5);	dhindt = BCs(6);
	velw	=	mwat/(rhow*AwaterX);

	double	P, L1, hout, mdotA, Tt1, Tt2, Tw1, Tw2;
	P	= State(1);		L1	=	State(2);	hout = State(4);	mdotA	=	State(5);
	Tt1	= State(7);		Tt2	=	State(8);
	Tw1 = State(10);	Tw2	=	State(11);

	double dPdt, dhoutdt, dL1dt, dTt1dt, dTt2dt, dTw1dt, dTw2dt, detA;

	if(!zone1.setstateP(P))
	{
		errorlog.Add("MBShellTubeEvaporator::fnTPSHDerivatives: ","Could not set saturated properties.");
		return false;
	}

	hv		=	zone1.gethv();
	rhov	=	zone1.getrhov();

	h1bar = (hrin + hv)/2.0;
	h2bar = (hv + hout)/2.0;

	if(!zone1.setstate(P,h1bar))
	{
		errorlog.Add("MBShellTubeEvaporator::fnTPSHDerivatives: ","Could not set state in zone 1.");
		return false;
	}
	if(!zone2.setstate(P,h2bar))
	{
		errorlog.Add("MBShellTubeEvaporator::fnTPSHDerivatives: ","Could not set state in zone 2.");
		return false;
	}

	rho1bar	=	zone1.getrho();
	rho2bar	=	zone2.getrho();

	T1bar = zone1.getT();
	T2bar = zone2.getT();

	drhodh1	=	zone1.getdrdh();	drhodP1	=	zone1.getdrdP();
	drhodh2	=	zone2.getdrdh();	drhodP2	=	zone2.getdrdP();
	dhvdP	=	zone1.getdhvdP();

	velr1 = (fabs(mrin + mdotA)/2.0)/(rho1bar*farea);
	velr2 = (fabs(mdotA + mrout)/2.0)/(rho2bar*farea);

	alphar1 =	rhtcoeff(P,h1bar,velr1,Tt1);
	alphar2	=	rhtcoeff(P,h2bar,velr2,Tt2);

	alphaw1	=	whtcoeff(Tw1,velw);
	alphaw2	=	whtcoeff(Tw2,velw);

	Qr1 =	alphar1*Ahtout*L1*(T1bar - Tt1);
	Qr2	=	alphar2*Ahtout*(L - L1)*(T2bar - Tt2);
	dQflux = fabs(Qr1)/Ahtout;
	Qr	=	Qr1 + Qr2;
	Qw1 =	alphaw1*Ahtin*L1*(Tt1 - Tw1);
	Qw2	=	alphaw2*Ahtin*(L - L1)*(Tt2 - Tw2);
	Qw	=	Qw1 + Qw2;

	a11	=	farea*L1*(0.5*drhodh1*dhvdP + drhodP1);
	a12	=	farea*(rho1bar - rhov);
	a21	=	farea*L1*(0.5*rho1bar*dhvdP + 0.5*h1bar*drhodh1*dhvdP + h1bar*drhodP1 - 1.0);
	a22	=	farea*(rho1bar*h1bar - rhov*hv);
	a31	=	farea*(L - L1)*(0.5*drhodh2*dhvdP + drhodP2);
	a32	=	-farea*(rho2bar - rhov);
	a33	=	0.5*farea*(L - L1)*drhodh2;
	a41	=	farea*(L - L1)*(0.5*rho2bar*dhvdP + 0.5*h2bar*drhodh2*dhvdP + h2bar*drhodP2 - 1.0);
	a42	=	-farea*(rho2bar*h2bar - rhov*hv);
	a43	=	0.5*farea*(L - L1)*(rho2bar + h2bar*drhodh2);

	b11	=	a11 + a31;		b12 =	a12 + a32;		b13 =	a33;
	b21	=	a21 + hv*a31;	b22	=	a22 + hv*a32;	b23	=	hv*a33;
	b31	=	a41 - hv*a31;	b32	=	a42 - hv*a32;	b33	=	a43 - hv*a33;

	d1	=	mrin - mrout - 0.5*farea*L1*drhodh1*dhindt;
	d2	=	mrin*hrin - mrout*hv - Qr1 - 0.5*farea*L1*(rho1bar + h1bar*drhodh1)*dhindt;
	d3	=	mrout*(hv - hout) - Qr2;

	detA=	b11*b22*b33-b11*b23*b32-b21*b12*b33+b21*b13*b32+b31*b12*b23-b31*b13*b22;

	dPdt	=	((b22*b33-b23*b32)*d1+(-b12*b33+b13*b32)*d2+(b12*b23-b13*b22)*d3)/detA;
	dL1dt	=	((-b21*b33+b23*b31)*d1+(b11*b33-b13*b31)*d2+(-b11*b23+b13*b21)*d3)/detA;
	dhoutdt	=	((b21*b32-b22*b31)*d1+(-b11*b32+b12*b31)*d2+(b11*b22-b12*b21)*d3)/detA;

	dTt1dt	=	(Qr1 - Qw1 - TubeCap*(Tt1 - Tt2))/(TubeCap*L1);
	dTt2dt	=	(Qr2 - Qw2)/(TubeCap*(L - L1));
	dTw1dt	=	(Qw1 + mwat*Cpw*(Tw2 - Tw1) - WatCap*(Tw1 - Tw2)*dL1dt)/(WatCap*L1);
	dTw2dt	=	(Qw2 + mwat*Cpw*(Twin - Tw2))/(WatCap*(L - L1));

	mdotA	=	mrout + a31*dPdt + a32*dL1dt + a33*dhoutdt;

	Derivs(1,dPdt);		Derivs(2,dL1dt);	Derivs(4,dhoutdt);
	Derivs(5,mdotA);	Derivs(7,dTt1dt);	Derivs(8,dTt2dt);
	Derivs(10,dTw1dt);	Derivs(11,dTw2dt);

	return true;
}//MBShellTubeEvaporator::fnTPSHDerivatives

bool	MBShellTubeEvaporator::fnDerivatives(Vector& bcs, Vector& state, Vector& derivs,double& Qr, double& Qw)
{
	switch(iOPMODE)
	{
	case	TP:
		if(!fnTPDerivatives(bcs,state,derivs,Qr,Qw))
		{
			errorlog.Add("MBShellTubeEvaporator::fnSH: ","Could not compute derivatives.");
			return false;
		}
		break;
	case	TPSH:
		if(!fnTPSHDerivatives(bcs,state,derivs,Qr,Qw))
		{
			errorlog.Add("MBShellTubeEvaporator::fnSHTP: ","Could not compute derivatives.");
			return false;
		}
		break;
	}//switch iOPMODE

	return true;
}//MBShellTubeEvaporator::fnDerivatives

bool	MBShellTubeEvaporator::fnIntegrateEXPE(Vector& bcs, Vector& state, double dt)
{
	Vector	derivs;
	derivs.make(12);
	double Qr, Qw;

	//State variables
	double	P, L1, L2, hout, mdotA, mdotB, Tt1, Tt2, Tt3, Tw1, Tw2, Tw3;
	P		= state(1);		L1		= state(2); 	L2		= state(3); 
	hout	= state(4);		mdotA	= state(5);		mdotB	= state(6);
	Tt1		= state(7);		Tt2		= state(8);		Tt3		= state(9);
	Tw1		= state(10);	Tw2		= state(11);	Tw3		= state(12);

	double dPdt, dL1dt, dL2dt, dhoutdt, dTt1dt, dTt2dt, dTt3dt, dTw1dt, dTw2dt, dTw3dt;

	//Compute derivatives
	if(!fnDerivatives(bcs,state,derivs,Qr,Qw))
	{
		errorlog.Add("MBShellTubeEvaporator::fnIntegrateEXPE: ","Could not compute derivatives.");
		return false;
	}
	dPdt	=	derivs(1);		dL1dt	=	derivs(2);		dL2dt	=	derivs(3);
	dhoutdt	=	derivs(4);		mdotA	=	derivs(5);		mdotB	=	derivs(6);
	dTt1dt	=	derivs(7);		dTt2dt	=	derivs(8);		dTt3dt	=	derivs(9);
	dTw1dt	=	derivs(10);		dTw2dt	=	derivs(11);		dTw3dt	=	derivs(12);

	P	+=	dPdt*dt;		L1	+=	dL1dt*dt;	L2	+=	dL2dt*dt;
	hout	+=	dhoutdt*dt;
	Tt1 +=	dTt1dt*dt;		Tt2	+=	dTt2dt*dt;	Tt3	+=	dTt3dt*dt;
	Tw1 +=	dTw1dt*dt;		Tw2	+=	dTw2dt*dt;	Tw3	+=	dTw3dt*dt;

	state(1,P);		state(2,L1);	state(3,L2);
	state(4,hout);	state(5,mdotA);	state(6,mdotB);
	state(7,Tt1);	state(8,Tt2);	state(9,Tt3);
	state(10,Tw1);	state(11,Tw2);	state(12,Tw3);

	double mrin, hrin, mrout, rmass, rengy;
	mrin = bcs(3);	mrout = bcs(5);
	hrin = bcs(4);
	fnComputeRefMassEngy(rmass,rengy);
	dDYNREFQTY	=	rmass;
	dTOTMASSIN	+=	mrin*dt;
	dTOTMASSOUT	+=	mrout*dt;
	dREFMASSIMB	=	(dDYNREFQTY - dORIGREFQTY)-(dTOTMASSIN-dTOTMASSOUT);

	dDYNREFENGY	=	rengy;
	dTOTENGYIN	+=	mrin*hrin*dt;
	dTOTENGYOUT	+=	mrout*hout*dt + Qr;
	dENGYIMB	=	(dDYNREFENGY - dORIGREFENGY) - (dTOTENGYIN - dTOTENGYOUT);

	return true;
}//MBShellTubeEvaporator::fnIntegrateEXPE

bool	MBShellTubeEvaporator::fnChkModeChange(Vector& NewState, int& NewMode)
{
	r134astate	exitstate;
	double P, hout, hv, hl;
	P	=	NewState(1);
	hout=	NewState(4);
	if(!exitstate.setstateP(P))
	{
		errorlog.Add("MBShellTubeEvaporator::fnChkModeChange: ","Could not set exit state.");
		return false;
	}
	hv = exitstate.gethv();		hl = exitstate.gethl();
	switch(iOPMODE)
	{
	case	TP:
		if(hout > hv)				//output is TP
		{			NewMode	=	TPSH;		iMODECHANGE = true;	iMODECHANGETYPE	=	TP_TPSH;		}
		else
		{			NewMode	=	iOPMODE;	iMODECHANGE = false;									}
		break;
	case	TPSH:
		if(hout < hv)				//output is SH
		{			NewMode	=	TP;			iMODECHANGE = true;	iMODECHANGETYPE	=	TPSH_TP;		}
		else
		{			NewMode	=	iOPMODE;	iMODECHANGE = false;									}
		break;
	}

	return true;
}//MBShellTubeEvaporator::fnChkModeChange

bool	MBShellTubeEvaporator::fnAdvance1step(double& dt, double dhindt)
{
	//Boundary conditions
	double	mwat, Twin, mrin, hrin, mrout;
	mwat = VINPUTS(1);	Twin = VINPUTS(2); mrin = VINPUTS(3); hrin = VINPUTS(4);
	mrout= VINPUTS(5);
	//State variables
	double	P, L1, L2, hout, mdotA, mdotB, Tt1, Tt2, Tt3, Tw1, Tw2, Tw3;
	P		= VState(1);	L1		= VState(2); 	L2		= VState(3); 
	hout	= VState(4);	mdotA	= VState(5);	mdotB	= VState(6);
	Tt1		= VState(7);	Tt2		= VState(8);	Tt3		= VState(9);
	Tw1		= VState(10);	Tw2		= VState(11);	Tw3		= VState(12);
	double	L;
	L = VGeometry(3);

	Vector bcs, state, derivs;
	bcs.make(6);	state.make(12);	derivs.make(12);
	bcs(1,mwat);	bcs(2,Twin);	bcs(3,mrin);	bcs(4,hrin);	bcs(5,mrout);	bcs(6,dhindt);
	state = VState;
	int	NewMode;
	double	epsl;
	epsl = EPSL0;

	if(!fnIntegrateEXPE(bcs,state,dt))
	{
		errorlog.Add("MBShellTubeEvaporator::fnAdvance1step: ","Could not integrate.");
		return false;
	}
	if(!fnChkModeChange(state,NewMode))
		//SPEED IMPROVEMENT IDEA...pass hv, hl from Integrate to ChkModeChange
		//to avoid re-evaluation of state properties
	{
		errorlog.Add("MBShellTubeEvaporator::fnAdvance1step: ","Could not check mode.");
		return false;
	}
	if(iMODECHANGE)
	{
		VState = state;
		iOPMODE	=	NewMode;
		if(iMODECHANGETYPE == TP_TPSH)
		{
			VState(2,L - epsl);	VState(5,mrout);
			VState(8,Tt1);	VState(11,Tw1);
		}
		iMODECHANGE = false;
	}//if iMODECHANGE
	else
		VState = state;
	VOUTPUTS = VState;

	return true;
}//MBShellTubeEvaporator::fnAdvance1step

bool	MBShellTubeEvaporator::fnAdvance1sec()
{
	char errstr[512];

	double tstart, tend, dt;
	tstart = 0.0;	tend = 1.0;
	dt = dTSTEP;

	r134astate	exitstate, instate;
	double	exitx, inx;

	double dhindt;
	dhindt = fnDhindt(VINPUTS(4));
	printf("dhindt = %f%c",dhindt,13);

	while(tstart < tend)
	{
		if(!fnAdvance1step(dt,dhindt))
		{
			sprintf(errstr,"Could not advance evaporator state at time %f",tstart);
			errorlog.Add("MBShellTubeEvaporator::fnAdvance1sec: ",errstr);
			return false;
		}
		instate.setstate(VState(1),VINPUTS(4));
		inx = instate.getx();
		exitstate.setstate(VState(1),VState(4));
		exitx = exitstate.getx();
		tstart	+=	dt;
		dTIME	+=	dt;
		dt = (tstart + dt)>tend?(tend - tstart):dTSTEP;
	}

	return true;
}//MBShellTubeEvaporator::fnAdvance1sec