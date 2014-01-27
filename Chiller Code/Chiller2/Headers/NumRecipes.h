#ifndef	_NUMRECIPES
#define	_NUMRECIPES
#include "Headers\Matrix.h"

void	nrerror(const char *error_text);
float	*vector(int nl,int nh);
int		*ivector(int nl, int nh);
double	*dvector(int nl, int nh);
float	**matrix(int nrl,int nrh,int ncl,int nch);
double	**dmatrix(int nrl, int nrh, int ncl, int nch);
int		**imatrix(int nrl, int nrh, int ncl, int nch);
double	**submatrix(double **a, int oldrl, int oldrh, int oldcl, int oldch, int newrl, int newcl);
void	free_vector(float *v,int nl,int nh);
void	free_vector(int *v,int nl,int nh);
void	free_vector(double *v, int nl, int nh);
void	free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void	free_matrix(double **m, int nrl, int nrh, int ncl, int nch);
void	free_matrix(int **m, int nrl, int nrh, int ncl, int nch);
void	free_submatrix(double **b, int nrl, int nrh, int ncl, int nch);
float	**convert_matrix(float *a, int nrl, int nrh, int ncl, int nch);
void	free_convert_matrix(float **b, int nrl, int nrh, int ncl, int nch);
void	spline(Vector& x,Vector& y,int n, double yp1,double ypn, Vector& y2);
double	splint(Vector& xa, Vector& ya, Vector& y2a, int n, double x);
void	splie2(Vector& x1a, Vector& x2a, Matrix& ya, int m, int n, Matrix& y2a);

double	splin2(Vector& x1a, Vector& x2a, Matrix& ya, Matrix& y2a, int m, int n, double x1, double x2);
void	ludcmp(Matrix& a, int n, Vector& indx, double *d);
void	lubksb(Matrix& a, int n, Vector& indx, Vector& b);
void	mprove(Matrix& a, Matrix& alud, int n, Vector& indx, Vector& b, Vector& x);
double	roundoff(double, int);
int		SuccOverRelax(Matrix& A, Vector& B, Vector& X, double tol, int iter, double omega, double& dtmax);
bool	GaussSeidel(Matrix& A, Vector& B, Vector& X);

bool	TDMA1D(Matrix& A, Vector& B, Vector& X);
bool	TDMA2D(Matrix& AP, Matrix& AN, Matrix& AS, Matrix& AE, Matrix& AW, Matrix& B, Matrix& X);
void	Pow(Matrix& A, int pwr, Matrix& PwrA);
void	Exp(Matrix& A, Matrix& ExpA);
void	Eye(Matrix& A);

bool	Invert(Matrix& A, Matrix& Ainv);
bool	solveLU(Matrix& A, Vector& B, Vector& X);
bool	solveLU(Matrix& A, Vector& B, Vector& X, int n);

#endif