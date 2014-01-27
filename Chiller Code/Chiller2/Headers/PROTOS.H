#ifndef		_PROTOS
#define		_PROTOS

#define		EVAP 1
#define		COND 2

#define		MYZERO	1e-10
#define		INF		1e30
#define		MYPI	3.1415926535897932384626433832795
#define		MINIMAL	0
#define		FULL	1
#define		SAVE	2
#define		RESTORE	3

//Controller response states
#define		WAIT	0
#define		READY	1

//Controller algorithm modes
#define		STARTUP		1	//
#define		NORMAL		2	//Modes of operation
#define		SHUTDOWN	3	//

//Refrigerant states
#define		SUBCOOLED	1	//
#define		TWOPHASE	2	//Refrigerant phases
#define		SUPERHEATED	3	//

//Tolerances during iterative convergence
#define		TTOL	1e-3	//C
#define		STOL	1e-6	//kJ/kg-K
#define		MTOL	0.025	//kg
#define		HTOL	1e-3	//kJ/kg
#define		FTOL	1e-2	//Function residual
#define		MAXIMB	0.025	//Maximum allowable refrigerant imbalance
#define		MINIMB	0.01	//Minimum imbalance before tstep grows

//Limits
#define		ETAPMIN	0.30	//Minimum polytropic efficiency
#define		DELTA	0.001	//Central differencing step-size
#define		DEFDT	0.010	//Default step size for the first SUTIME seconds
#define		TCSTEP	0.040	//Default integration step-size for condenser
#define		TESTEP	0.010	//Default integration step-size for evaporator
#define		MAXDT	0.2		//Maximum integration step-size
#define		MINDT	0.002	//Minimum integration step-size
#define		DTINC	0.005	//tstep growth
#define		MAXITERS 5000	//Maximum number of iterations
#define		SUTIME	150		//Startup time s

//Refrigerants
#define		R134A	1
#define		R12		2
#define		R22		3

//Integration algorithms
#define		EXPE	1
#define		EXPEPC	2
#define		RK4		3
#define		IMPEPC	4
#define		IMPE	5		//Equivalent to the sequential solution method

#define		BULBLOCN	12	//Number of the evaporator node in which the bulb is located
#define		HINDEPTH	2	//Number of seconds to look back in computing dhindt

#define		EPSL0		0.01		//Least value of L1 or L2 at mode change during early start-up
#define		EPSL1		0.003	//Least value of L1 or L2 at mode change during normal operation

#define		SH		1		//All superheated
#define		TP		2		//All two-phase
#define		SHTP	3		//Superheated at entry, two-phase at exit
#define		TPSH	4		//Two-phase at entry, superheated at exit
#define		TPSC	5		//Two-phase at entry, sub-cooled at exit
#define		SCTP	6		//Sub-cooled at entry, two-phase at exit
#define		SHTPSC	7		//Superheated at entry, two-phase intermediate, sub-cooled at exit
#define		SCTPSH	8		//Sub-cooled at entry, two-phase intermediate, superheated at exit

#define		SH_SHTP			1
#define		SHTP_SHTPSC		2
#define		SHTPSC_SHTP		3
#define		SHTP_SH			4
#define		TP_TPSH			5
#define		TPSH_TP			6

#endif