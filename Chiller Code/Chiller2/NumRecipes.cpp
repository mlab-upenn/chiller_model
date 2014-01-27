#include	"StdAfx.h"
#include	<stdlib.h>
#include	<string.h>
#include	<math.h>
#include	"Headers\Matrix.h"
#include	"Headers\NumRecipes.h"
#include "Headers\ErrorLog.h"

#define		MYZERO	1e-10
#define		INF		1e30
#define		MYPI	3.1415926535897932384626433832795
#define		LOG2INV	1.4426950408889634073599246810019

extern ErrorLog	errorlog;

float *vector(int nl,int nh)
//Allocates a float vector with range [nl...nh]
{
	float *v;

	v	=	new float[nh-nl+1];
	// (float *)malloc((unsigned)(nh-nl+1)*sizeof(float));
	if(!v)
	{
		errorlog.Add("vector() ","Memory allocation failed");
		return NULL;
	}
	return v-nl;
}

int *ivector(int nl, int nh)
//Allocates an integer vector with range [nl...nh]
{
	int *v;

	v	=	new int[nh-nl+1];
	//malloc((unsigned)(nh-nl+1)*sizeof(int));
	if(!v)
	{
		errorlog.Add("ivector() ","Memory allocation failed");
		return NULL;
	}
	return v-nl;
}

double *dvector(int nl, int nh)
//Allocates a double vector with range [nl...nh]
{
	double *v;

	v	=	new double[nh-nl+1];
	//malloc((unsigned)(nh-nl+1)*sizeof(double));
	if(!v)
	{
		errorlog.Add("dvector() ","Memory allocation failed");
		return NULL;
	}
	return v-nl;
}

float **matrix(int nrl,int nrh,int ncl,int nch)
//Allocates a float matrix with range [nrl...nrh][ncl...nch]
{
	int i;
	float **m;

	//Allocate pointers to rows
	m	=	new float*[nrh-nrl+1];
	//*)[malloc((unsigned)(nrh-nrl+1)*sizeof(float*));
	if(!m)
	{
		errorlog.Add("matrix() ","Memory allocation failed");
		return NULL;
	}
	m -= nrl;

	//Allocate rows and set pointers to them
	for(i=nrl;i<=nrh;i++)
	{
		m[i]	=	new float[nch-ncl+1];
		//*)malloc((unsigned)(nch-ncl+1)*sizeof(float));
		if(!m[i])
		{
			errorlog.Add("matrix() ","Memory allocation failed");
			return NULL;
		}
		m[i] -= ncl;
	}

	//Return pointer to array of pointers to rows
	return m;
}

double **dmatrix(int nrl, int nrh, int ncl, int nch)
//Allocates a double matrix with range [nrl...nrh][ncl...nch]
{
	int i;
	double **m;

	//Allocate pointers to rows
	m	=	new double*[nrh-nrl+1];
	//*)malloc((unsigned)(nrh-nrl+1)*sizeof(double*));
	if(!m)
	{
		errorlog.Add("dmatrix() ","Memory allocation failed");
		return NULL;
	}
	m -= nrl;

	//Allocate rows and set pointers to them
	for(i=nrl;i<=nrh;i++)
	{
		m[i]	=	new double[nch-ncl+1];
		//*)malloc((unsigned)(nch-ncl+1)*sizeof(double));
		if(!m[i])
		{
			errorlog.Add("dmatrix() ","Memory allocation failed");
			return NULL;
		}
		m[i] -= ncl;
	}
	//Return pointer to array of pointers to rows
	return m;
}

int **imatrix(int nrl, int nrh, int ncl, int nch)
//Allocates an int matrix with range [nrl...nrh][ncl...nch]
{
	int i, **m;

	//Allocate pointers to rows
	m	=	new int*[nrh-nrl+1];
	//*malloc((unsigned)(nrh-nrl+1)*sizeof(int*));
	if(!m)
	{
		errorlog.Add("imatrix() ","Memory allocation failed");
		return NULL;
	}
	m -= nrl;

	//Allocate rows and set pointers to them
	for(i=nrl;i<=nrh;i++)
	{
		m[i]	=	new int[nch-ncl+1];
		//*) malloc((unsigned)(nch-ncl+1)*sizeof(int*));
		if(!m[i])
		{
			errorlog.Add("imatrix() ","Memory allocation failed");
			return NULL;
		}
		m[i] -= ncl;
	}
	//Return pointer to array of pointers to rows
	return m;
}

double **submatrix(double **a, int oldrl, int oldrh, int oldcl, int oldch, int newrl, int newcl)
//Returns a submatrix with range [newrl...newrl+(oldrh-oldrl)][newcl...newcl+(oldch-oldcl)]
//pointing to the existing matrix range a[oldrl...oldrh][oldcl...oldch]
{
	int i,j;
	double **m;

	//Allocate pointers to rows
	m	=	new double*[oldrh-oldrl+1];
	//*)malloc((unsigned)(oldrh-oldrl+1)*sizeof(double*));
	if(!m)
	{
		errorlog.Add("submatrix() ","Memory allocation failed");
		return NULL;
	}
	m -= newrl;

	//Set pointers to rows
	for(i=oldrl,j=newrl; i<=oldrh;i++,j++)
		m[j] = a[i]+oldcl-newcl;

	//Return pointer to array of pointers to rows
	return m;
}

void free_vector(float *v,int nl,int nh)
//Frees a float vector allocated by vector()
{
	float *vptr = v+nl;
	if(vptr==NULL) 
		return;
	delete []vptr;
	vptr=NULL;
}

void free_vector(int *v,int nl,int nh)
//Frees a int vector allocated by ivector()
{
	int *vptr = v+nl;
	if(vptr==NULL)
		return;
	delete []vptr;
	vptr=NULL;
}

void free_vector(double *v, int nl, int nh)
//Frees a double vector allocated by dvector()
{
	double *vptr = v+nl;
	if(vptr==NULL)
		return;
	delete []vptr;
	vptr=NULL;
}

void free_matrix(float **m, int nrl, int nrh, int ncl, int nch)
//Frees a matrix allocated with matrix()
{
	int i;
	float **mptr=m+nrl, *m1ptr;

	if(mptr==NULL)
		return;
	for(i=nrh;i>=nrl;i--)
	{
		m1ptr = m[i]+ncl;
		if(m1ptr==NULL)
			return;
		delete []m1ptr;
		m1ptr = NULL;
	}
	delete []mptr;
	mptr = NULL;
}

void free_matrix(double **m, int nrl, int nrh, int ncl, int nch)
//Frees a matrix allocated with dmatrix()
{
	int i;
	double **mptr=m+nrl, *m1ptr;

	if(mptr==NULL)
		return;
	for(i=nrh;i>=nrl;i--)
	{
		m1ptr = m[i]+ncl;
		if(m1ptr==NULL)
			return;
		delete []m1ptr;
		m1ptr = NULL;
	}
	delete []mptr;
	mptr = NULL;
}

void free_matrix(int **m, int nrl, int nrh, int ncl, int nch)
//Frees a matrix allocated with imatrix()
{
	int i;
	int **mptr=m+nrl, *m1ptr;

	if(mptr==NULL)
		return;
	for(i=nrh;i>=nrl;i--)
	{
		m1ptr = m[i]+ncl;
		if(m1ptr==NULL)
			return;
		delete []m1ptr;
		m1ptr = NULL;
	}
	delete []mptr;
	mptr = NULL;
}

void free_submatrix(double **b, int nrl, int nrh, int ncl, int nch)
//Frees a submatrix allocated by submatrix()
{
	double **bptr = b+nrl;
	if(bptr==NULL)
		return;
	delete [](b+nrl);
	bptr = NULL;
}

float **convert_matrix(float *a, int nrl, int nrh, int ncl, int nch)
//Allocate a float matrix m[nrl...nrh][ncl...nch] that points to the matrix a declared int the standard C manner
// as a[nrow][ncol], where nrow = nrh-nrl+1 and ncol = nch-ncl+1.  The routine should be called with the address
//&a[0][0] as the first argument
{
	int i,j,nrow,ncol;
	float **m;

	nrow = nrh-nrl+1;
	ncol = nch-ncl+1;

	//Allocate pointers to rows
	m = new float*[nrow];
	//*)malloc((unsigned)(nrow)*sizeof(float*));
	if(!m)
	{
		errorlog.Add("convert_matrix() ","Memory allocation failed");
		return NULL;
	}
	m -= nrl;
	for(i=0,j=nrl;i<=nrow-1;i++,j++)
		m[j] = a+ncol*i-ncl;	//Set pointers to rows
	return m;					//Return pointer to array of pointers to rows
}

void free_convert_matrix(float **b, int nrl, int nrh, int ncl, int nch)
//Frees a matrix allocated by convert_matrix()
{
	float **bptr = b+nrl;
	if(bptr==NULL)
		return;
	delete [](b+nrl);
	bptr = NULL;
}

void spline(Vector& x,Vector& y,int n, double yp1,double ypn, Vector& y2)
//Given arrays x[1...n] and y[1...n] containing a tabulated function, i.e, y_i = f(x_i), with x_1<x_2<...x_n
// and given alues yp1 and ypn for the first derivative of the interpolating function at points 1 and n
//respectively, this routine returns an array y2[1...n] that contains the second derivatives of the interpolating
//function at the tabulated points x_i.  If yp1 and/or ypn are equal to 1e30 or larger, the routine is signalled
//to set the corresponding boundary condition for a natural spline, with zero second derivative on that boundary
{
	int i,k;
	double p,qn,sig,un;
	Vector u;

	u.make(n-1);
	if(!u.ok())
		return;
	if(yp1>0.99e30)	//The lower boundary condition is set either to be "natural"
	{
		y2(1,0.0);
		u(1,0.0);
	}
	else			//or else to have a specified first derivative
	{
		y2(1,-0.5);
		u(1,(3.0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1));
	}
	for(i=2;i<=n-1;i++)
	{
		sig = (x(i)-x(i-1))/(x(i+1)-x(i-1));
		p	= sig*y2(i-1)+2.0;
		y2(i,(sig-1.0)/p);
		u(i,(y(i+1)-y(i))/(x(i+1)-x(i)) - (y(i)-y(i-1))/(x(i)-x(i-1)));
		u(i,(6.0*u(i)/(x(i+1)-x(i-1))-sig*u(i-1))/p);
	}
	if(ypn>0.99e30)	//The upper boundary condition is set either to be "natural"
		qn = un = 0.0;
	else			//or else to have a specified first derivative
	{
		qn = 0.5;
		un = (3.0/(x(n)-x(n-1)))*(ypn - (y(n)-y(n-1))/(x(n)-x(n-1)));
	}
	y2(n,(un-qn*u(n-1))/(qn*y2(n-1)+1.0));
	for(k=n-1;k>=1;k--)
		y2(k,y2(k)*y2(k+1)+u(k));
}

double splint(Vector& xa, Vector& ya, Vector& y2a, int n, double x)
//Given the arrays xa[1...n] and ya[1...n], which tabulate a function (with the xa_i's in order), and given the
//array y2a[1...n], which is the output from 'spline', and given a value of x, this routine returns a cubic-
//spline interpolated value y
{
	int klo, khi, k;
	double h,b,a;

	klo = 1;
	khi = n;

	while(khi-klo>1)
	{
		k = (khi+klo) >> 1;
		if(xa(k)>x)
			khi = k;
		else
			klo = k;
	}
	h = xa(khi)-xa(klo);
	if(h==0.0)
	{
		errorlog.Add("splint ","Independant variable column contains duplicate entries");
		return -1;
	}
	if(x>xa(n)||x<xa(1))
		errorlog.AddWarning("splint ","","Extrapolation detected");
	a = (xa(khi)-x)/h;
	b = (x-xa(klo))/h;

	return a*ya(klo) + b*ya(khi)+((a*a*a-a)*y2a(klo)+(b*b*b-b)*y2a(khi))*(h*h)/6.0;
}

void splie2(Vector& x1a, Vector& x2a, Matrix& ya, int m, int n, Matrix& y2a)
//Given a tabulated function ya[1...m][1...n] and tabulated independent variables x1a[1...m] and x2a[1...n],
//this routine constructs one-dimensional natural cubic splines of the rows of ya and returns the second
//derivatives in the array y2a[1...m][1...n]
{
	int j;
	Vector yaj, y2aj(m);

	for(j=1;j<=m;j++)
	{
		ya.getrow(j,yaj);
		spline(x2a,yaj,n,1e30,1e30,y2aj);
		y2a.putrow(j,y2aj);
	}
}

double splin2(Vector& x1a, Vector& x2a, Matrix& ya, Matrix& y2a, int m, int n, double x1, double x2)
//Given x1a, x2a, ya, m,n, as described in splie2, and y2a as produced by that routine, and, given a desired
//interpolating point x1, x2, this routine returns an interpolated function value y by bicubic spline interpolation
{
	int j;
	Vector ytmp, yytmp, y2aj, yaj;

	ytmp.make(m);
	yytmp.make(m);
	y2aj.make(n);
	yaj.make(n);
	if(!ytmp.ok()||!yytmp.ok())
		return -1;

	for(j=1;j<=m;j++)
	{
		ya.getrow(j,yaj);
		y2a.getrow(j,y2aj);
		yytmp(j,splint(x2a,yaj,y2aj,n,x2));
	}
	spline(x1a,yytmp,m,1e30,1e30,ytmp);
	return splint(x1a,yytmp,ytmp,m,x1);
}

void ludcmp(Matrix& a, int n, Vector& indx, double *d)
//Given an n x n matrix a[1...n][1...n], this routine replaces it by the LU decomposition of a rowwise permutation 
//of itself.  a and n are input.  a is output.  indx[1...n] is an output vector which records the row permutation 
//effected by the partial pivoting.  d is output as +-1, depending on whether the number of row interchanges was 
//even or odd, respectively.  This routine is used in combination with lubksb to solve linear equations or invert
//a matrix
{
	int i, imax, j,k;
	double big, dum, sum, temp;
	Vector vv;		//vv stores the implicit scaling of each row

	vv.make(n);
	if(!vv.ok())
		return;
	*d = 1.0;			//no row interchanges yet
	for(i=1;i<=n;i++)	//Loop over rows to get the implicit scaling information
	{
		big = 0.0;
		for(j=1;j<=n;j++)
			if((temp=fabs(a(i,j)))>big)
				big = temp;
		if(big==0.0)
		{
			errorlog.Add("ludcmp ","Matrix is singular");
			return;
		}
		vv(i,1.0/big);	//save the scaling
	}
	for(j=1;j<=n;j++)		//This is the loop over columns of Crout's method
	{
		for(i=1;i<j;i++)
		{
			sum = a(i,j);
			for(k=1;k<i;k++)
				sum -= a(i,k)*a(k,j);
			a(i,j,sum);
		}
		big = 0.0;		//initialize for the search for the largest pivot element
		for(i=j;i<=n;i++)
		{
			sum = a(i,j);
			for(k=1;k<j;k++)
				sum -= a(i,k)*a(k,j);
			a(i,j,sum);
			if((dum=vv(i)*fabs(sum))>=big)	//is the figure of merit for the pivot
			{								//better than the best so far?
				big = dum;
				imax = i;
			}
		}
		if(j!=imax)				//If we need to interchange rows...
		{
			for(k=1;k<=n;k++)
			{
				dum = a(imax,k);
				a(imax,k,a(j,k));
				a(j,k,dum);
			}
			*d = -(*d);			//...and change the parity of d
			vv(imax,vv(j));	//also interchange the scale factor
		}
		indx(j,double(imax));
		if(a(j,j)==0.0)
			a(j,j,MYZERO);
		//If the pivot element is zero the matrix is singular (at least to the precision of the algorithm)
		//For some applications on singular matrices, it is desirable to substitute MYZERO for zero
		if(j!=n)		//Now, finally, divide by the pivot element
		{
			dum = 1.0/(a(j,j));
			for(i=j+1;i<=n;i++)
				a(i,j,a(i,j)*dum);
		}
	}	//Go back for the next column in the reduction
}//ludcmp

void lubksb(Matrix&  a, int n, Vector& indx, Vector& b)
//Solves the set of n linear equations AX=B.  Here, a[1...n][1...n] is input, not as the matrix A, but as its
//LU decomposition, determined by the routine ludcmp.  indx[1...n] is input as the permutation vector returned
//by ludcmp.  b[1...n] is input as the right-hand-side vector B, and returns with the solution vector X.  a,n 
//and indx are not modified by this routine and can be left in place for successive calls with different right
//hand sides b.  This routine takes into account the possibility that b will begin with many zero elements, so
//it is efficient for use in matrix inversion.
{
	int i, ii=0,ip,j;
	double sum;

	for(i=1;i<=n;i++)
	{
		ip = (int)indx(i);
		sum = b(ip);
		b(ip,b(i));
		if(ii)
			for(j=ii;j<=i-1;j++)
				sum -= a(i,j)*b(j);
		else if(sum)
			ii = i;
		b(i,sum);
	}
	for(i=n;i>=1;i--)
	{
		sum = b(i);
		for(j=i+1;j<=n;j++)
			sum -= a(i,j)*b(j);
		b(i,sum/a(i,i));
	}
}

void	mprove(Matrix& a, Matrix& alud, int n, Vector& indx, Vector& b, Vector& x)
{
	int i,j;
	double sdp;
	Vector r;
	r.make(n);
	for(i=1;i<=n;i++)
	{
		sdp	=	-b(i);
		for(j=1;j<=n;j++)
			sdp	+=	a(i,j)*x(j);
		r(i,sdp);
	}
	lubksb(alud,n,indx,r);
	for(i=1;i<=n;i++)
		x(i,x(i)-r(i));
}

double	roundoff(double num, int dec)
{
	//Rounds off num to dec decimal places
	double t1, t2, thalf, mult=1.0;

	mult = pow(10,dec);
	t1	=	floor(num*mult);
	t2	=	ceil(num*mult);
	thalf = (t1+t2)/2.0;
	if((num*mult)>thalf)
		return t2/mult;
	else
		return t1/mult;
}

bool	TDMA1D(Matrix& A, Vector& B, Vector& X)
{
	//One dimensional, Tri-Diagonal Matrix Algorithm
	//...A is the coefficient matrix, of size Nx3 where N is the number of cells
	//...B is the source vector
	//...X is the state vector
	//As passed, the discrete equations are required to be in the form
	// a_p.phi_P = a_w*phi_W + a_e*phi_E + B_p

	//Check input data for dimensional consistency
	bool GoodArgs = false;
	GoodArgs = (A.cols()==3)&&(A.rows()==B.size())&&(B.size()==X.size());
	if(!GoodArgs)
	{
		errorlog.Add("TDMA1D ","Dimensions of input arguments are inconsistent.");
		return false;
	}
	//Data format checked and found OK

	//Check A for singularity...
	for(int i=1;i<=A.rows();i++)
	{
		if(fabs(A(i,2))<=MYZERO)
		{
			errorlog.Add("TDMA1D ","Zero entry found on main diagonal of A.");
			return false;
		}
	}//'A' Matrix is healthy

	//Proceed with solution...
	int N;
	double multiplier;
	N = A.rows();
	//...fix A to suit required format of equation structure viz. AX = B
	for(i=1;i<=N;i++)
	{
		A(i,1,-1.0*A(i,1));
		A(i,3,-1.0*A(i,3));
	}

	//Forward elimination...
	for(i=2;i<=N;i++)
	{
		multiplier	= A(i,1)/A(i-1,2);
		A(i,1,multiplier);
		A(i,2,A(i,2) - multiplier*A(i-1,3));
		B(i,B(i) - multiplier*B(i-1));
	}

	//Back substitution...
	X(N,B(N)/A(N,2));
	for(i=N-1;i>=1;i--)
		X(i,(B(i) - A(i,3)*X(i+1))/A(i,2));

	//Undo changes made to A
	for(i=1;i<=N;i++)
	{
		A(i,1,-1.0*A(i,1));
		A(i,3,-1.0*A(i,3));
	}

	return true;
}

bool	TDMA2D(Matrix& AP, Matrix& AN, Matrix& AS, Matrix& AE, Matrix& AW, Matrix& B, Matrix& X)
{
	//Two dimensional, Tri-Diagonal Matrix Algorithm
	//...AP, AN, AS, AE and AW are the coefficient matrices, all are NIxNJ sized matrices
	//...X is the NIxNJ sized state matrix, 
	//... passed in as a guess and output as the refined X after one to-and-fro sweep
	//....each in vertical and horizontal directions

	//As passed, the discrete equations are required to be in the form
	// a_p.phi_P = a_w*phi_W + a_e*phi_E + a_n*phi_N + a_s*phi_S + B_p


	//Check input data for dimensional consistency
	bool GoodArgs = false;
	GoodArgs = (AP.rows()==AN.rows())&&(AN.rows()==AS.rows())&&(AS.rows()==AE.rows())&&(AE.rows()==AW.rows());
	GoodArgs = GoodArgs&&(AP.cols()==AN.cols())&&(AN.cols()==AS.cols())&&(AS.cols()==AE.cols())&&(AE.cols()==AW.cols());
	GoodArgs = GoodArgs&&(AP.rows()==X.rows())&&(AP.cols()==X.cols());
	if(!GoodArgs)
	{
		errorlog.Add("TDMA2D ","Inconsistent dimensions found in input matrices.");
		return false;
	}//Input arguments are of consistent dimensions

	//Proceed with solution...
	int NI, NJ;
	NI = AP.cols();		//Get number of horizontal discretizations
	NJ = AP.rows();		//Get number of vertical discretizations

	Matrix Aj;		Aj.make(NI,3);
	Vector Bj, Xj;	Bj.make(NI);	Xj.make(NI);

	//Begin sweeps...
	//Begin upward vertical sweep - through j = constant lines from j = 1 to NJ
	double bji;
	for(int j=1;j<=NJ;j++)
	{
		Aj.reset();	Bj.reset();	Xj.reset(); bji = 0.0;
		//Construct A and B structures for TDMA1D
		for(int i=1;i<=NI;i++)
		{
			Aj(i,1,AW(j,i));
			Aj(i,2,AP(j,i));
			Aj(i,3,AE(j,i));
			bji = B(j,i);
			if(j==1)
				bji	+= AN(j,i)*X(j+1,i);
			else if(j==NJ)
				bji += AS(j,i)*X(j-1,i);
			else
				bji += AN(j,i)*X(j+1,i) + AS(j,i)*X(j-1,i);
			Bj(i,bji);
		}//for i
		//Solve for jth row using TDMA1D
		if(!TDMA1D(Aj,Bj,Xj))
		{
			errorlog.Add("TDMA2D:TDMA1D ","TDMA1D Failed.  See error log file for more information.");
			return false;
		}
		//Update state matrix with most recently computed states
		X.putrow(j,Xj);
	}//for j - End upward vertical sweep
	//Begin downward vertical sweep - through j = constant lines from j = NJ to 1
	for(j=NJ;j>=1;j--)
	{
		Aj.reset();	Bj.reset();	Xj.reset(); bji = 0.0;
		//Construct A and B structures for TDMA1D
		for(int i=1;i<=NI;i++)
		{
			Aj(i,1,AW(j,i));
			Aj(i,2,AP(j,i));
			Aj(i,3,AE(j,i));
			bji = B(j,i);
			if(j==1)
				bji	+= AN(j,i)*X(j+1,i);
			else if(j==NJ)
				bji += AS(j,i)*X(j-1,i);
			else
				bji += AN(j,i)*X(j+1,i) + AS(j,i)*X(j-1,i);
			Bj(i,bji);
		}//for i
		//Solve for jth row using TDMA1D
		if(!TDMA1D(Aj,Bj,Xj))
		{
			errorlog.Add("TDMA2D:TDMA1D ","TDMA1D Failed.  See error log file for more information.");
			return false;
		}
		//Update state matrix with most recently computed states
		X.putrow(j,Xj);
	}//for j - End downward vertical sweep
	Aj.make(NJ,3);
	Bj.make(NJ);	Xj.make(NJ);
	//Begin right horizontal sweep - through i = constant lines from i = 1 to NI
	for(int i=1;i<=NI;i++)
	{
		Aj.reset();	Bj.reset();	Xj.reset(); bji = 0.0;
		//Construct A and B structures for TDMA1D
		for(j=1;j<=NJ;j++)
		{
			Aj(j,1,AS(j,i));
			Aj(j,2,AP(j,i));
			Aj(j,3,AN(j,i));
			bji = B(j,i);
			if(i==1)
				bji	+= AE(j,i)*X(j,i+1);
			else if(i==NI)
				bji += AW(j,i)*X(j,i-1);
			else
				bji += AE(j,i)*X(j,i+1) + AW(j,i)*X(j,i-1);
			Bj(j,bji);
		}//for i
		//Solve for jth row using TDMA1D
		if(!TDMA1D(Aj,Bj,Xj))
		{
			errorlog.Add("TDMA2D:TDMA1D ","TDMA1D Failed.  See error log file for more information.");
			return false;
		}
		//Update state vector with most recently computed states
		X.putcol(i,Xj);
	}//End right horizontal sweep
	//Begin left horizontal sweep - through i = constant lines from i = NI to 1
	for(i=NI;i>=1;i--)
	{
		Aj.reset();	Bj.reset();	Xj.reset(); bji = 0.0;
		//Construct A and B structures for TDMA1D
		for(j=1;j<=NJ;j++)
		{
			Aj(j,1,AS(j,i));
			Aj(j,2,AP(j,i));
			Aj(j,3,AN(j,i));
			bji = B(j,i);
			if(i==1)
				bji	+= AE(j,i)*X(j,i+1);
			else if(i==NI)
				bji += AW(j,i)*X(j,i-1);
			else
				bji += AE(j,i)*X(j,i+1) + AW(j,i)*X(j,i-1);
			Bj(j,bji);
		}//for i
		//Solve for jth row using TDMA1D
		if(!TDMA1D(Aj,Bj,Xj))
		{
			errorlog.Add("TDMA2D:TDMA1D ","TDMA1D Failed.  See error log file for more information.");
			return false;
		}
		//Update state vector with most recently computed states
		X.putcol(i,Xj);
	}//End left horizontal sweep

	return true;
}

int SuccOverRelax(Matrix& A, Vector& B, Vector& X, double tol, int iter, double omega, double& dtmax)
{
	Vector R, deltaX;
	double	maxdeltaX=0.0;
	int i,j, n,iterno=0,maxiters;
	bool DONE = false;

	n	=	A.rows();
	maxiters = iter;
	R.make(n);
	deltaX.make(n);

	while(!DONE&&iterno<=maxiters)
	{
		iterno++;
		for(i=1;i<=n;i++)
		{
			for(j=1;j<=n;j++)
				R(i,R(i)+A(i,j)*X(j));
			R(i,B(i) - R(i));
			deltaX(i,omega*R(i)/A(i,i));
			X(i,X(i) + deltaX(i));
			R(i,0.0);
		}
		maxdeltaX = 0.0;
		for(i=1;i<=n;i++)
			maxdeltaX = (fabs(deltaX(i))>maxdeltaX)?fabs(deltaX(i)):maxdeltaX;
		if(maxdeltaX<tol)
		{
			DONE = true;
			dtmax = maxdeltaX;
		}
	}//while

	return iterno;
}

bool	GaussSeidel(Matrix& A, Vector& B, Vector& X)
{
	//Gauss Siedel iterative solution of a system of linear algebraic equations
	//represented by AX = B
	//A -> System matrix
	//B -> Source vector
	//X -> State vector...contains initial guess values on entry
	//		and converged values on exit

	//Check arguments for size and format
	bool GoodArgs = false;
	GoodArgs = (A.rows()==A.cols())&&(A.rows()==B.size())&&(B.size()==X.size());
	if(!GoodArgs)
	{
		errorlog.Add("GaussSeidel ","Badly formatted arguments passed to Gauss Siedel routine.");
		return false;
	}

	//Ensure non-zero diagonal entries
	for(int i=1;i<=A.rows();i++)
	{
		if(fabs(A(i,i))<=MYZERO)
		{
			errorlog.Add("GaussSeidel ","Zero entry found on main diagonal of A.");
			return false;
		}
	}//'A' Matrix is healthy

	//Proceed with solution...
	int j, n, iters=0;
	double res, Change, OldNorm;
	n = A.rows();	//Determine number of cells

	//...the iteration
	do
	{
		OldNorm = X.norm();		//Record norm of X
		for(i=1;i<=n;i++)
		{
			res = 0.0;
			for(j=1;j<=n;j++)
				if(j!=i)	res += A(i,j)*X(j);
			X(i,(B(i)-res)/A(i,i));
		}
		Change = fabs(X.norm() - OldNorm);	//Compute change in X
		iters++;
		//Check for failure to converge
		if(iters>1000)
		{
			errorlog.Add("GaussSeidel ","Convergence not achieved in 1000 iterations.");
			return false;
		}
	}while(Change>1e-6); //while Change is high

	return true;
}

bool solveLU(Matrix& A, Vector& B, Vector& X)
{
	double d;
	int n;
	Vector indx;
	Matrix a, alud;
	Vector b, blud;

	n = A.rows();

	a.make(n,n);	alud.make(n,n);
	b.make(n);		blud.make(n);
	indx.make(n);
	
	if(errorlog.IsError())	
	{	
		errorlog.Add("solveLU ","Temporary variables could not be created");	
		return false;	
	}

	a = A;	alud = A;
	b = B;	blud = B;

	ludcmp(alud,n,indx,&d);
	lubksb(alud,n,indx,blud);
	mprove(a,alud,n,indx,b,blud);
	mprove(a,alud,n,indx,b,blud);

	X = blud;
	return true;
}

bool solveLU(Matrix& A, Vector& B, Vector& X, int n)
{
	double d;
	Vector indx;
	Matrix a, alud;
	Vector b, blud;

	a.make(n,n);	alud.make(n,n);
	b.make(n);		blud.make(n);
	indx.make(n);
	
	if(errorlog.IsError())	
	{	
		errorlog.Add("solveLU ","Temporary variables could not be created");	
		return false;	
	}

	for(int i=1;i<=n;i++)
	{
		for(int j=1;j<=n;j++)
		{
			a(i,j,A(i,j));
			alud(i,j,A(i,j));
		}
		b(i,B(i));
		blud(i,B(i));
	}

	ludcmp(alud,n,indx,&d);
	lubksb(alud,n,indx,blud);
	mprove(a,alud,n,indx,b,blud);
	mprove(a,alud,n,indx,b,blud);

	for(i=1;i<=n;i++)
		X(i,blud(i));
	return true;
}

bool	Invert(Matrix& A, Matrix& Ainv)
{
	int	n;
	n = A.rows();
	Vector indx, cols;
	indx.make(n);
	cols.make(n);
	double	d;
	Matrix Alud;	Alud.make(n,n);
	Alud = A;

	ludcmp(Alud, n, indx, &d);
	for(int j=1;j<=n;j++)
	{
		cols.reset();
		cols(j,1.0);
		lubksb(Alud,n,indx,cols);
		for(int i=1;i<=n;i++)
			Ainv(i,j,cols(i));
	}

	return true;
}

void	Eye(Matrix& A)
{
	int M;
	M = A.rows();
	A.reset();
	for(int i=1;i<=M;i++)
		A(i,i,1.0);
}

void	Pow(Matrix& A, int pwr, Matrix& PwrA)
{
	PwrA = A;
	for(int i=2;i<=pwr;i++)
		PwrA = PwrA*A;
}

void	Exp(Matrix&A, Matrix& ExpA)
{
	double k;
	k = log(A.RowInfinityNorm())*LOG2INV;
	k = int(k)+1;
	double scale;
	scale = pow(2,k);
	A = A/scale;

	int L;
	double facti=1.0;
	Matrix Ap;
	Ap.make(A.rows(),A.cols());
	Matrix ExpAp;
	ExpAp.make(A.rows(),A.cols());
	L = A.ExpSeriesLen();
	Eye(ExpAp);
	ExpAp = ExpAp + A;

	for(int i=2;i<=L;i++)
	{
		facti *= double(i);
		Pow(A,i,Ap);
		ExpAp = ExpAp + Ap/facti;
	}
	for(i=1;i<=k;i++)
	{
		Pow(ExpAp,2,ExpA);
		ExpAp = ExpA;
	}

	A = A*scale;
}
