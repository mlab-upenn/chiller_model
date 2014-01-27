#ifndef _COMPONENTS
#define	_COMPONENTS

#include	"Headers\Matrix.h"

class	Component
{
public:
	Vector				VGeometry;
	Vector*				MGeometry;
	Vector				VState;
	Vector*				MState;
	virtual void		fnSetState		()		const = 0;
	virtual bool		fnReEvalState	()		const = 0;
	Component							()		{}
	virtual bool		fnDefineGeometry()		const = 0;
	virtual bool		fnAdvance1sec	()		const = 0;
};
//************************************************************
class	FVShellTubeHX : public Component
{
private:
	int		iHXTYPE, iGCOLS,iSCOLS;
	double	dORIGREFQTY, dDYNREFQTY, dTOTMASSIN, dTOTMASSOUT, dREFMASSIMB;
	double	dTSTEP;
	double	TotVolume;
	bool	bTSIZING;

	virtual void	fnSetState			()		const	{}
	virtual bool	fnDefineGeometry	()		const	{return true;}
	virtual bool	fnReEvalState		()		const	{return true;}
	virtual bool	fnAdvance1sec		()		const	{return true;}
	bool			fnReEvalStateExpE	(double&,Matrix&,double&);
	bool			fnReEvalStateExpEPC	(double&,Matrix&,double&);
	bool			fnReEvalStateRK4	(double&,Matrix&,double&);
	bool			fnReEvalStateImpEPC	(double&,Matrix&,double&);
	bool			fnReEvalStateImpE	(double&,Matrix&,double&);

	//Sub-functions used exclusively in the sequential solution
	double			fnRefMassBalance(double mr0g, double rhog, double rho, double tstep);
	double			fnRefEnergyBalance(double hgi_1, double hgi1, double Ttgi, double Trgi, double mgj, 
									   double mgj_1, double rho, double hi, double Hr, double tstep);
	double			fnTubeEnergyBalance(double Hr, double Trg, double Hw, double Twg, double Tt,double tstep);
	double			fnWaterEnergyBalance(double Twgi1, double Ttgi, double Tw, double Hw, double tstep);
	bool			fnNodalMassEnergySolver(int nodeno, double Pg, double rho, Matrix& newstate2, double tstep);
	bool			fnHXMassEnergyBalancesSweep(double Pg, Vector& rho, Matrix& newstate2, double tstep);
	bool			fnEnthalpyProfileIterator(double Pg, Vector& rho, Matrix& newstate2, double tstep);
	bool			fnPressureIterator(double& newstate1, Matrix& newstate2, double& HXEnBal, double tstep);
	///////////////

	double			fnRefHTCoeff		(Vector& proplist,double Tt);
	double			fnWatHTCoeff		(Vector&);
	bool			FVShellTubeHX::fnPopulateAB(double mrin, double hrin, double mrout,Vector& hi, Vector& Qri);

public:
	int		iN_INPUTS, iN_OUTPUTS;
	Vector	VINPUTS, VOUTPUTS;
	int		NNODES;
	int		NTUBES;
	double	dTIME;
	Matrix	A;
	Vector	ais, bis, cis, dis;
	Vector	partsumsa;
	Vector	rpropertylist;
	Vector	wpropertylist;
	Vector	X, B;
	double	sumQ, sumc;

	FVShellTubeHX();
	int		NODES();
	double	Volume();
	bool	fnDefineGeometry		(const char*);
	bool	fnSetState				(double, double, Matrix&);
	bool	fnSetState				(double, double, double);
	void	fnGetState				(double&,Matrix&);
	void	fnGetState				(Vector&);
	bool	fnStateDerivs			(double,Vector&,Vector&,Vector&,Vector&,Vector&,double&,Vector&,Vector&,Vector&,double&,double&);
	bool	fnStateDerivs			(double,Vector&,Vector&,Vector&,double&,Vector&,Vector&,Vector&);
	bool	fnAdvance1sec			(int);
};

//************************************************************

class	MBShellTubeCondenser : public Component
{
private:
	int		iHXTYPE, iGCOLS,iSCOLS;
	double	dORIGREFQTY, dDYNREFQTY, dTOTMASSIN, dTOTMASSOUT, dREFMASSIMB;
	double	dORIGREFENGY,dDYNREFENGY,dTOTENGYIN, dTOTENGYOUT, dENGYIMB;
	double	dTSTEP;
	double	TotVolume;
	bool	bTSIZING;
	int		iOPMODE;
	int		iMODECHANGE;
	int		iMODECHANGETYPE;
	Vector	dHRINOLD;

	virtual void	fnSetState			()		const	{}
	virtual bool	fnDefineGeometry	()		const	{return true;}
	virtual bool	fnReEvalState		()		const	{return true;}
	virtual bool	fnAdvance1sec		()		const	{return true;}
	bool			fnReEvalState		(double&,Matrix&,double&);

	double			rhtcoeff			(double P, double h, double vel, double Tt);
	double			whtcoeff			(double Tw, double vel);
	double			fnRefHTCoeff		(Vector& proplist,double Tt);
	double			fnWatHTCoeff		(Vector&);
	bool			fnMeanPropsSC		(double P, double h1, double h2, double dL, double& rhobar, double& rhoubar);
	bool			fnMeanPropsSH		(double P, double h1, double h2, double dL, double& rhobar, double& rhoubar);
	bool			fnMeanPropsTP		(double P, double x1, double x2, double& rhobar, double& rhoubar);
	void			fnComputeRefMassEngy(double& mass, double& engy);
	double			fnDhindt			(double hrin);

public:
	int		iN_INPUTS, iN_OUTPUTS;
	Vector	VINPUTS, VOUTPUTS;
	int		NTUBES;
	double	dTIME;
	Vector	rpropertylist;
	Vector	wpropertylist;

	MBShellTubeCondenser();
	bool	fnDefineGeometry		(const char*);
	bool	fnInitialize			(Vector&);

	bool	fnSH					();
	bool	fnTP					();
	bool	fnSHTP					();
	bool	fnTPSC					();
	bool	fnSHTPSC				();
	bool	fnSHDerivatives			(Vector&,Vector&,Vector&,double& Qr, double& Qw);
	bool	fnTPDerivatives			(Vector&,Vector&,Vector&,double& Qr, double& Qw);
	bool	fnSHTPDerivatives		(Vector&,Vector&,Vector&,double& Qr, double& Qw);
	bool	fnTPSCDerivatives		(Vector&,Vector&,Vector&,double& Qr, double& Qw);
	bool	fnSHTPSCDerivatives		(Vector&,Vector&,Vector&,double& Qr, double& Qw);
	bool	fnDerivatives			(Vector&,Vector&,Vector&,double& Qr, double& Qw);
	bool	fnIntegrateEXPE			(Vector&,Vector&,double);
	bool	fnIntegrateEPC			(Vector&,Vector&,double);
	bool	fnAdvance1step			(double&,double);
	bool	fnComputeDt				(double,double,double,double,double&);
	bool	fnAdvance1sec			();
	bool	fnChkModeChange			(Vector&,int&);

	int		fnOpMode				()	{return iOPMODE;}
	double	fnRefMass				()	{return dDYNREFQTY;}
	double	fnRefMassImbalance		()	{return dREFMASSIMB;}

};

//************************************************************

class	MBShellTubeEvaporator : public Component
{
private:
	int		iHXTYPE, iGCOLS,iSCOLS;
	double	dORIGREFQTY, dDYNREFQTY, dTOTMASSIN, dTOTMASSOUT, dREFMASSIMB;
	double	dORIGREFENGY,dDYNREFENGY,dTOTENGYIN, dTOTENGYOUT, dENGYIMB;
	double	dTSTEP;
	double	TotVolume;
	bool	bTSIZING;
	int		iOPMODE;
	int		iMODECHANGE;
	int		iMODECHANGETYPE;
	Vector	dHRINOLD;
	double	dQflux;

	virtual void	fnSetState			()		const	{}
	virtual bool	fnDefineGeometry	()		const	{return true;}
	virtual bool	fnReEvalState		()		const	{return true;}
	virtual bool	fnAdvance1sec		()		const	{return true;}
	bool			fnReEvalState		(double&,Matrix&,double&);

	double			rhtcoeff			(double P, double h, double vel, double Tt);
	double			whtcoeff			(double Tw, double vel);
	double			fnRefHTCoeff		(Vector& proplist,double Tt);
	double			fnWatHTCoeff		(Vector&);
	bool			fnMeanPropsSH		(double P, double h1, double h2, double dL, double& rhobar, double& rhoubar);
	bool			fnMeanPropsTP		(double P, double x1, double x2, double& rhobar, double& rhoubar);
	void			fnComputeRefMassEngy(double& mass, double& engy);
	double			fnDhindt			(double hrin);

public:
	int		iN_INPUTS, iN_OUTPUTS;
	Vector	VINPUTS, VOUTPUTS;
	int		NTUBES;
	double	dTIME;
	Vector	rpropertylist;
	Vector	wpropertylist;

	MBShellTubeEvaporator();
	bool	fnDefineGeometry		(const char*);
	bool	fnInitialize			(Vector&);

	bool	fnTP					();
	bool	fnTPSH					();
	bool	fnTPDerivatives			(Vector&,Vector&,Vector&,double& Qr, double& Qw);
	bool	fnTPSHDerivatives		(Vector&,Vector&,Vector&,double& Qr, double& Qw);
	bool	fnDerivatives			(Vector&,Vector&,Vector&,double& Qr, double& Qw);
	bool	fnIntegrateEXPE			(Vector&, Vector&, double);
	bool	fnAdvance1step			(double&,double);;
	bool	fnComputeDt				(double,double,double,double,double&);
	bool	fnAdvance1sec			();
	bool	fnChkModeChange			(Vector&,int&);

	int		fnOpMode				()	{return iOPMODE;}
	double	fnRefMass				()	{return dDYNREFQTY;}
	double	fnRefMassImbalance		()	{return dREFMASSIMB;}
};
//************************************************************

class	MicroTech	:	public Component
{
private:
	int		iRUNMODE;
	int		iACTIONMODE;
	double	dACTIONTIME, dTIMER;
	bool	bOPEN, bCLOSE;
	double	dDNSTEPFACTOR, dUPSTEPFACTOR;
	virtual	void	fnSetState			()		const	{}
	virtual bool	fnDefineGeometry	()		const	{return true;}
	virtual bool	fnReEvalState		()		const	{return true;}
	virtual bool	fnAdvance1sec		()		const	{return true;}
	double			fnReEvalState		(Vector&);

public:
	MicroTech();
	void	fnDefineGeometry			(Vector&);
	void	fnSetState					(double,int);
	void	fnGetState					(double&);
	bool	fnAdvance1sec				(Vector&,double&);
};

//************************************************************
class	CentComp : public Component
{
private:
	int iMODE;
	virtual	void	fnSetState			()		const	{}
	virtual bool	fnDefineGeometry	()		const	{return true;}
	virtual bool	fnReEvalState		()		const	{return true;}
	virtual bool	fnAdvance1sec		()		const	{return true;}
	bool			fnReEvalState		(Vector&);
	MicroTech		Controller;

	double	fnWOVFlowRate				(double,double,double);
	double	fnPolytropicEfficiency		(double,double);
	double	fnRatedLoadAmpsPerCent		(double);
	double	fnInvRatedLoadAmpsPerCent	(double);
	double	fnPolytropicWork			(double,double,double,double);

	bool	fnQuasiSteadyState1D		(double,double,double,double,double,double&);
	bool	fnQuasiSteadyState2D		(double P1, double h1, double P2, double Pwr, double h2g, double mdotrg, double& h2, double& mdotr);
	bool	fnResidue1					(double,double,double,double,double,double&);
	bool	fnResidue2					(double,double,double,double,double&);

public:
	double	dTIME;
	int	iN_INPUTS,	iN_OUTPUTS;	
	Vector	VINPUTS, VOUTPUTS;
	CentComp();
	bool	fnDefineGeometry			(const char*);
	bool	fnSetState					(double,int,Vector&);
	bool	fnSetState					(double);
	void	fnGetState					(int&,Vector&);
	bool	fnAdvance1sec				();
	void	shutoff();
};
//************************************************************

class	LumpedCapacitance	: public Component
{
private:
	virtual	void	fnSetState			()		const	{}
	virtual bool	fnDefineGeometry	()	const	{return true;}
	virtual bool	fnReEvalState		()		const	{return true;}
	virtual bool	fnAdvance1sec		()		const	{return true;}
	double			fnReEvalState		(double);
public:
	LumpedCapacitance();
	bool	fnDefineGeometry	(double);
	void	fnSetState			(double);
	void	fnGetState			(double&);
	bool	fnAdvance1sec		(double,double&);
};

//************************************************************
class	ThermoStaticExValve : public Component
{
private:
	virtual	void	fnSetState		()		const	{}
	virtual bool	fnDefineGeometry()		const	{return true;}
	virtual bool	fnReEvalState	()		const	{return true;}
	virtual bool	fnAdvance1sec	()		const	{return true;}
	bool			fnReEvalState1	();		//Original physical model
	bool			fnReEvalState2	();		//Map only model
	bool			fnReEvalState3	();		//Updated model with bulb-temperature map
	bool			fnReEvalState4	();		//Updated model with bulb-temperature map
	bool			fnReEvalState5	(double);		//Updated model with bulb-temperature map

public:
	double	dTIME;
	int	iN_INPUTS, iN_OUTPUTS;
	Vector	VINPUTS, VOUTPUTS;

	ThermoStaticExValve();
	bool			fnDefineGeometry	(const char*);
	bool			fnSetState			(double,double,double);
	bool			fnSetState			(double);
	void			fnGetState			(Vector&);
	virtual bool	fnAdvance1sec		(double);
};

//************************************************************

class	Orifice	:	public Component
{
private:
	virtual	void	fnSetState		()		const	{}
	virtual bool	fnDefineGeometry()		const	{return true;}
	virtual bool	fnReEvalState	()		const	{return true;}
	virtual bool	fnAdvance1sec	()		const	{return true;}
	bool			fnReEvalState	(double&);
public:
	int iN_INPUTS, iN_OUTPUTS;
	Vector	VINPUTS, VOUTPUTS;
	Orifice();
	bool	fnDefineGeometry		(const char*);
	bool	fnDefineGeometry		(double);
	void	fnSetState				(double);
	void	fnGetState				(double&);
	bool	fnAdvance1sec			();
};

//************************************************************

class	VapCompCentLiqChiller	:	public Component
{
private:
	virtual	void	fnSetState		()		const	{}
	virtual bool	fnDefineGeometry()		const	{return true;}
	virtual bool	fnReEvalState	()		const	{return true;}
	virtual bool	fnAdvance1sec	()		const	{return true;}
	bool			fnReEvalState	(Vector&,Vector&);
	bool	bINITIALIZED;

	FVShellTubeHX*			Evaporator;
	FVShellTubeHX*			Condenser;
	CentComp*				Compressor;
	ThermoStaticExValve*	Valve;
	Orifice*				CoolingLineBypass;

	double	dNCMass;

public:

	double	dTIME;
	bool	bNONCONDENSABLES;
	int		iN_INPUTS, iN_OUTPUTS;
	int		iNNODES_C,	iNNODES_E;
	Vector	VINPUTS, VOUTPUTS;

	VapCompCentLiqChiller();
	bool	fnSetState			(double time, double devapP, Matrix& MevapS, 
								 double dcondP, Matrix& McondS, 
								 int imode, Vector& VcompS, 
								 double dmdotv, double dTbulb, double dmdotcl,double dTewo,double NCMass);
	bool	fnSetState			(double Tewo, double Tcwo, double Mtotal, double NCMass);
	bool	fnDefineGeometry	(const char*);
	bool	fnDefineGeometry	();
	bool	fnSaveState			(const char*);
	bool	fnSaveState			();
	bool	fnGetState			(Matrix&);
	bool	fnLoadState			(const char*);
	bool	fnLoadState			(int);
	bool	fnLoadState			();
	bool	fnAdvance1sec		(int);
	bool	fnInitialized		();
	void	shutoff				();

	void	fnSetNCMass			(double);
	double	fnGetNCMass			();

};


//************************************************************
#endif