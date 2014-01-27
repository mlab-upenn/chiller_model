#ifndef	_PROPERTIES
#define _PROPERTIES
#include "Headers\Matrix.h"

class two_phase_table
{
	Matrix	table;
	Matrix	splines;
	Vector	currprops;
	bool	PROPSRGOOD;

	int rows;
	int cols;

public:
	two_phase_table();
	bool load2ptab(char *infile);
	Vector& getprop(double x, Vector& props);
	bool	goodprops()	{	return PROPSRGOOD;	}
};

class one_phase_table
{
	Matrix	table;
	Matrix	splines;
	Vector	currprops;
	bool	PROPSRGOOD;

	int rows;
	int cols;
	int blkind[100];
	int nblks;
	int maxblk();

public:
	one_phase_table();
	bool load1ptab(char *infile);
	Vector& getprop(double x1, double x2, Vector& props);
	bool	goodprops()	{	return PROPSRGOOD;	}
};

class watstate
{
protected:
	double	T,		//Temperature C
			k,		//Conductivity, W/m.C
			rho,	//Density, kg/m^3
			mu,		//Viscosity, kg/m.s
			Cp;		//Specific heat, kJ/kg.C
	bool	UPDATE;
public:
	watstate();
	~watstate();
	bool	setstate(double Temp);
	double	getT();
	double	getk();
	double	getrho();
	double	getmu();
	double	getCp();
	bool	getupdate();
};

class r134astate
{
protected:
	double	P,		//Pressure	kPa
			T,		//Temperature C
			Ts,		//Saturation temperature C
			hl,		//Saturated liquid enthalpy		kJ/kg
			hv,		//Saturated vapor enthalpy		kJ/kg
			h,		//enthalpy						kJ/kg
			vl,		//saturated liquid specific volume	m^3/kg
			vv,		//saturated vapor specific volume	m^3/kg
			v,		//specific volume					m^3/kg
			rhol,	//saturated liquid density			kg/m^3
			rhov,	//saturated vapor density			kg/m^3
			rho,	//density							kg/m^3
			sl,		//saturated liquid entropy			kJ/kg.C
			sv,		//saturated vapor entropy			kJ/kg.C
			s,		//entropy							kJ/kg.C
			ul,		//saturated liquid internal energy	kJ/kg
			uv,		//saturated vapor internal energy	kJ/kg
			u,		//internal energy					kJ/kg
			kl,		//saturated liquid thermal conductivity	W/m.C
			kv,		//saturated vapor thermal conductivity	W/m.C
			k,		//thermal conductivity					W/m.C
			mul,	//saturated liquid viscosity			kg/m.s
			muv,	//saturaved vapor viscosity				kg/m.s
			mu,		//viscosity								kg/m.s
			Cpl,	//saturated liquid specific heat		kJ/kg.C
			Cpv,	//saturated vapor specific heat			kJ/kg.C
			Cp,		//specific heat							kJ/kg.C
			x,		//quality
			drdP,	//partial rho wrt P at constant h
			drdh,	//partial rho wrt h at constant P
			drvdP,	//derivative of rhov wrt P
			drldP,	//derivative of rhol wrt P
			dhvdP,	//derivative of hv wrt P
			dhldP;	//derivative of hl wrt P
	int		phase;	//phase	= 1 for subcooled, = 2 for two-phase and = 3 for superheated
	bool	UPDATE;	//flag to ensure properties are updated by a call to setstate
public:
	int		setphase(double enth, double hl, double hv);
	r134astate();
	~r134astate();

	bool	setstateP(double Press);
	bool	setstateT(double T);
	bool	setstateHS(double h, double s);
	bool	setstate(double P, double h);
	bool	setstatePV(double P, double v);
	bool	setstatePT(double P, double T);
	void	setP(double Pres)	{	P=Pres; 	UPDATE = false;	}
	void	seth(double enth)	{	h=enth; 	UPDATE = false;	}
	double	getP()			{	return P;		}
	double	getT()			{	return T;		}
	double	getTs()			{	return Ts;		}
	double	gethl()			{	return hl;		}
	double	gethv()			{	return hv;		}
	double	geth()			{	return h;		}
	double	getvl()			{	return vl;		}
	double	getvv()			{	return vv;		}
	double	getv()			{	return v;		}
	double	getrhol()		{	return rhol;	}
	double	getrhov()		{	return rhov;	}
	double	getrho()		{	return rho;		}
	double	getsl()			{	return sl;		}
	double	getsv()			{	return sv;		}
	double	gets()			{	return s;		}
	double	getul()			{	return ul;		}
	double	getuv()			{	return uv;		}
	double	getu()			{	return u;		}
	double	getkl()			{	return kl;		}
	double	getkv()			{	return kv;		}
	double	getk()			{	return k;		}
	double	getmul()		{	return mul;		}
	double	getmuv()		{	return muv;		}
	double	getmu()			{	return mu;		}
	double	getCpl()		{	return Cpl;		}
	double	getCpv()		{	return Cpv;		}
	double	getCp()			{	return Cp;		}
	double	getx()			{	return x;		}
	double	getdrdP()		{	return drdP;	}
	double	getdrdh()		{	return drdh;	}
	double	getdrvdP()		{	return drvdP;	}
	double	getdrldP()		{	return drldP;	}
	double	getdhvdP()		{	return dhvdP;	}
	double	getdhldP()		{	return dhldP;	}
	int		getphase()		{	return phase;	}
	bool	getupdate()		{	return UPDATE;	}
};//class r134astate

class r12state
{
protected:
	double	P,		//Pressure	kPa
			T,		//Temperature C
			Ts,		//Saturation temperature C
			hl,		//Saturated liquid enthalpy		kJ/kg
			hv,		//Saturated vapor enthalpy		kJ/kg
			h,		//enthalpy						kJ/kg
			vl,		//saturated liquid specific volume	m^3/kg
			vv,		//saturated vapor specific volume	m^3/kg
			v,		//specific volume					m^3/kg
			rhol,	//saturated liquid density			kg/m^3
			rhov,	//saturated vapor density			kg/m^3
			rho,	//density							kg/m^3
			sl,		//saturated liquid entropy			kJ/kg.C
			sv,		//saturated vapor entropy			kJ/kg.C
			s,		//entropy							kJ/kg.C
			ul,		//saturated liquid internal energy	kJ/kg
			uv,		//saturated vapor internal energy	kJ/kg
			u,		//internal energy					kJ/kg
			kl,		//saturated liquid thermal conductivity	W/m.C
			kv,		//saturated vapor thermal conductivity	W/m.C
			k,		//thermal conductivity					W/m.C
			mul,	//saturated liquid viscosity			kg/m.s
			muv,	//saturaved vapor viscosity				kg/m.s
			mu,		//viscosity								kg/m.s
			Cpl,	//saturated liquid specific heat		kJ/kg.C
			Cpv,	//saturated vapor specific heat			kJ/kg.C
			Cp,		//specific heat							kJ/kg.C
			x,		//quality
			drdP,	//partial rho wrt P at constant h
			drdh;	//partial rho wrt h at constant P
	int		phase;	//phase	= 1 for subcooled, = 2 for two-phase and = 3 for superheated
	bool	UPDATE;	//flag to ensure properties are updated by a call to setstate
public:
	int		setphase(double enth, double hl, double hv);
	r12state();
	~r12state();

	bool	setstateP(double Press);
	bool	setstateT(double T);
	bool	setstateHS(double h, double s);
	bool	setstate(double P, double h);
	bool	setstatePV(double P, double v);
	bool	setstatePT(double P, double T);
	void	setP(double Pres)	{	P=Pres; 	UPDATE = false;	}
	void	seth(double enth)	{	h=enth; 	UPDATE = false;	}
	double	getP()			{	return P;		}
	double	getT()			{	return T;		}
	double	getTs()			{	return Ts;		}
	double	gethl()			{	return hl;		}
	double	gethv()			{	return hv;		}
	double	geth()			{	return h;		}
	double	getvl()			{	return vl;		}
	double	getvv()			{	return vv;		}
	double	getv()			{	return v;		}
	double	getrhol()		{	return rhol;	}
	double	getrhov()		{	return rhov;	}
	double	getrho()		{	return rho;		}
	double	getsl()			{	return sl;		}
	double	getsv()			{	return sv;		}
	double	gets()			{	return s;		}
	double	getul()			{	return ul;		}
	double	getuv()			{	return uv;		}
	double	getu()			{	return u;		}
	double	getkl()			{	return kl;		}
	double	getkv()			{	return kv;		}
	double	getk()			{	return k;		}
	double	getmul()		{	return mul;		}
	double	getmuv()		{	return muv;		}
	double	getmu()			{	return mu;		}
	double	getCpl()		{	return Cpl;		}
	double	getCpv()		{	return Cpv;		}
	double	getCp()			{	return Cp;		}
	double	getx()			{	return x;		}
	double	getdrdP()		{	return drdP;	}
	double	getdrdh()		{	return drdh;	}
	int		getphase()		{	return phase;	}
	bool	getupdate()		{	return UPDATE;	}
};//class r12state

class r22state
{
protected:
	double	P,		//Pressure	kPa
			T,		//Temperature C
			Ts,		//Saturation temperature C
			hl,		//Saturated liquid enthalpy		kJ/kg
			hv,		//Saturated vapor enthalpy		kJ/kg
			h,		//enthalpy						kJ/kg
			vl,		//saturated liquid specific volume	m^3/kg
			vv,		//saturated vapor specific volume	m^3/kg
			v,		//specific volume					m^3/kg
			rhol,	//saturated liquid density			kg/m^3
			rhov,	//saturated vapor density			kg/m^3
			rho,	//density							kg/m^3
			sl,		//saturated liquid entropy			kJ/kg.C
			sv,		//saturated vapor entropy			kJ/kg.C
			s,		//entropy							kJ/kg.C
			ul,		//saturated liquid internal energy	kJ/kg
			uv,		//saturated vapor internal energy	kJ/kg
			u,		//internal energy					kJ/kg
			kl,		//saturated liquid thermal conductivity	W/m.C
			kv,		//saturated vapor thermal conductivity	W/m.C
			k,		//thermal conductivity					W/m.C
			mul,	//saturated liquid viscosity			kg/m.s
			muv,	//saturaved vapor viscosity				kg/m.s
			mu,		//viscosity								kg/m.s
			Cpl,	//saturated liquid specific heat		kJ/kg.C
			Cpv,	//saturated vapor specific heat			kJ/kg.C
			Cp,		//specific heat							kJ/kg.C
			x;		//quality
	int		phase;	//phase	= 1 for subcooled, = 2 for two-phase and = 3 for superheated
	bool	UPDATE;	//flag to ensure properties are updated by a call to setstate
public:
	int		setphase(double enth, double hl, double hv);
	r22state();
	~r22state();

	bool	setstateP(double Press);
	bool	setstateT(double T);
	bool	setstateHS(double h, double s);
	bool	setstate(double P, double h);
	bool	setstatePV(double P, double v);
	bool	setstatePT(double P, double T);
	void	setP(double Pres)	{	P=Pres; 	UPDATE = false;	}
	void	seth(double enth)	{	h=enth; 	UPDATE = false;	}
	double	getP()			{	return P;		}
	double	getT()			{	return T;		}
	double	getTs()			{	return Ts;		}
	double	gethl()			{	return hl;		}
	double	gethv()			{	return hv;		}
	double	geth()			{	return h;		}
	double	getvl()			{	return vl;		}
	double	getvv()			{	return vv;		}
	double	getv()			{	return v;		}
	double	getrhol()		{	return rhol;	}
	double	getrhov()		{	return rhov;	}
	double	getrho()		{	return rho;		}
	double	getsl()			{	return sl;		}
	double	getsv()			{	return sv;		}
	double	gets()			{	return s;		}
	double	getul()			{	return ul;		}
	double	getuv()			{	return uv;		}
	double	getu()			{	return u;		}
	double	getkl()			{	return kl;		}
	double	getkv()			{	return kv;		}
	double	getk()			{	return k;		}
	double	getmul()		{	return mul;		}
	double	getmuv()		{	return muv;		}
	double	getmu()			{	return mu;		}
	double	getCpl()		{	return Cpl;		}
	double	getCpv()		{	return Cpv;		}
	double	getCp()			{	return Cp;		}
	double	getx()			{	return x;		}
	int		getphase()		{	return phase;	}
	bool	getupdate()		{	return UPDATE;	}
};//class r22state

class r500state
{
protected:
	double	P,		//Pressure	kPa
			T,		//Temperature C
			Ts,		//Saturation temperature C
			hl,		//Saturated liquid enthalpy		kJ/kg
			hv,		//Saturated vapor enthalpy		kJ/kg
			h,		//enthalpy						kJ/kg
			vl,		//saturated liquid specific volume	m^3/kg
			vv,		//saturated vapor specific volume	m^3/kg
			v,		//specific volume					m^3/kg
			rhol,	//saturated liquid density			kg/m^3
			rhov,	//saturated vapor density			kg/m^3
			rho,	//density							kg/m^3
			sl,		//saturated liquid entropy			kJ/kg.C
			sv,		//saturated vapor entropy			kJ/kg.C
			s,		//entropy							kJ/kg.C
			ul,		//saturated liquid internal energy	kJ/kg
			uv,		//saturated vapor internal energy	kJ/kg
			u,		//internal energy					kJ/kg
			kl,		//saturated liquid thermal conductivity	W/m.C
			kv,		//saturated vapor thermal conductivity	W/m.C
			k,		//thermal conductivity					W/m.C
			mul,	//saturated liquid viscosity			kg/m.s
			muv,	//saturaved vapor viscosity				kg/m.s
			mu,		//viscosity								kg/m.s
			Cpl,	//saturated liquid specific heat		kJ/kg.C
			Cpv,	//saturated vapor specific heat			kJ/kg.C
			Cp,		//specific heat							kJ/kg.C
			x,		//quality
			drdP,	//partial rho wrt P at constant h
			drdh;	//partial rho wrt h at constant P
	int		phase;	//phase	= 1 for subcooled, = 2 for two-phase and = 3 for superheated
	bool	UPDATE;	//flag to ensure properties are updated by a call to setstate
	int		setphase(double enth, double hl, double hv);
public:
	r500state();
	~r500state();

	bool	setstateP(double Press);
	bool	setstateT(double T);
	void	setP(double Pres){	P=Pres; 	UPDATE = false;	}
	void	seth(double enth){	h=enth; 	UPDATE = false;	}
	double	getP()			{	return P;	}
	double	getT()			{	return T;	}
	double	getTs()			{	return Ts;	}
	double	gethl()			{	return hl;	}
	double	gethv()			{	return hv;	}
	double	geth()			{	return h;	}
	double	getvl()			{	return vl;	}
	double	getvv()			{	return vv;	}
	double	getv()			{	return v;	}
	double	getrhol()		{	return rhol;}
	double	getrhov()		{	return rhov;}
	double	getrho()		{	return rho;	}
	double	getsl()			{	return sl;	}
	double	getsv()			{	return sv;	}
	double	gets()			{	return s;	}
	double	getul()			{	return ul;	}
	double	getuv()			{	return uv;	}
	double	getu()			{	return u;	}
	double	getkl()			{	return kl;	}
	double	getkv()			{	return kv;	}
	double	getk()			{	return k;	}
	double	getmul()		{	return mul;	}
	double	getmuv()		{	return muv;	}
	double	getmu()			{	return mu;	}
	double	getCpl()		{	return Cpl;	}
	double	getCpv()		{	return Cpv;	}
	double	getCp()			{	return Cp;	}
	double	getx()			{	return x;	}
	double	getdrdP()		{	return drdP;}
	double	getdrdh()		{	return drdh;}
	int		getphase()		{	return phase;	}
	bool	getupdate()		{	return UPDATE;	}
};//class r500state

bool loadtables(int,bool);
void unloadtables(void);

#endif
