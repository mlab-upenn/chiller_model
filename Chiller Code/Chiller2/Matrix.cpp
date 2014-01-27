#include	"StdAfx.h"
#include	<iostream>
#include	<string>
#include	<math.h>
#include	"Headers\NumRecipes.h"
#include	"Headers\Matrix.h"
#include	"Headers\ErrorLog.h"

extern ErrorLog	errorlog;

Vector::Vector()
{
	if(vec!=NULL)
		vec = NULL;
	OK	= false;
}

Vector::Vector(int a)
{
	N = a;

	vec = dvector(1,a);
	if(vec==NULL)
	{
		OK = false;
		errorlog.Add("Vector::Vector ","Construction failed");
		return;
	}
	for(int i=1;i<=N;i++)
		vec[i] = 0.0;
	OK = true;
}
Vector::~Vector()
{
	if(vec==NULL)
		return;
	free_vector(vec,1,N);
}

bool	Vector::make(int a)
{
	if(OK&&vec!=NULL)	//if vector is already defined, redefine
		free_vector(vec,1,N);
	N = a;
	vec = dvector(1,a);
	if(vec==NULL)
	{
		OK = false;
		errorlog.Add("Vector::make ","Memory allocation failed");
		return OK;
	}
	for(int i=1;i<=N;i++)
		vec[i] = 0.0;
	OK = true;
	return OK;
}

int		Vector::size()						{	return	N;}
bool	Vector::ok()						{	return OK;}

Vector&	Vector::operator = (Vector& V)
{
	if(N!=V.size())
		return *this;
	for(int i=1;i<=N;i++)
		vec[i] = V(i);
	return *this;
}

Vector& Vector::operator + (Vector& V)
{
	if(N!=V.size())
		return *this;
	for(int i=1;i<=N;i++)
		vec[i] = vec[i] + V(i);
	return *this;
}

Vector& Vector::operator - (Vector& V)
{
	if(N!=V.size())
		return *this;
	for(int i=1;i<=N;i++)
		vec[i] = vec[i] - V(i);
	return *this;
}

void	Vector::reset()
{
	for(int i=1;i<=N;i++)
		vec[i] = 0.0;
}

double	Vector::min()
{
	double least;
	least = vec[1];
	for(int i=2;i<=N;i++)
	{
		if(vec[i]<least)
			least = vec[i];
	}
	return least;
}

double	Vector::max()
{
	double most;
	most = vec[1];
	for(int i=2;i<=N;i++)
	{
		if(vec[i]>most)
			most = vec[i];
	}
	return most;
}

double	Vector::sum()
{
	double total=0.0;
	int i=1;
	while(i<=N)
		total += vec[i++];
	return total;
}

double	Vector::norm()
{
	double total=0.0;
	int i=1;
	while(i<=N)
		total += vec[i]*vec[i++];
	return total;
}

double	Vector::fmin()
{
	double	least;
	least	=	fabs(vec[1]);
	for(int i=2;i<=N;i++)
	{
		if(fabs(vec[i])<least)
			least = fabs(vec[i]);
	}
	return least;
}

double	Vector::fmax()
{
	double most;
	most = fabs(vec[1]);
	for(int i=2;i<=N;i++)
	{
		if(fabs(vec[i])>most)
			most = fabs(vec[i]);
	}
	return most;
}

void	Vector::scale(double factor)
{
	for(int i=1;i<=N;i++)
		vec[i] = factor*vec[i];
}

//===============================MATRIX CLASS BEGINS======================================

Matrix::Matrix()	//default constructor - initializes mat to a 3x3 zero matrix
{
	N	=	0;		//initialize to a 3x3 matrix
	M	=	0;
	mat	=	NULL;
	DEF	=	false;
}//Matrix

Matrix::Matrix(int P)	//alternate constructor - initializes mat to a PxP zero matrix
{
	int i,j;
	N	=	P;
	M	=	P;
	mat	=	dmatrix(1,M,1,N);
	if(mat==NULL)
	{
		DEF = false;
		return;
	}
	DEF	=	true;
	for(i=1;i<=M;i++)
		for(j=1;j<=N;j++)
			mat[i][j] = 0.0;
}//Matrix(P)

Matrix::Matrix(int P,int Q)		//alternate constructor - initializes mat to a PxQ zero matrix
{
	int i,j;
	N	=	Q;
	M	=	P;
	mat	=	dmatrix(1,M,1,N);
	if(mat==NULL)
	{
		DEF = false;
		errorlog.Add("Matrix::Matrix ","Construction failed");
		return;
	}
	DEF	=	true;
	for(i=1;i<=M;i++)
		for(j=1;j<=N;j++)
			mat[i][j] = 0.0;
}//Matrix(P,Q)

bool	Matrix::make(int P, int Q)
{
	if(DEF)		//if defined, redefine
		free_matrix(mat,1,M,1,N);

	int i,j;
	N	=	Q;
	M	=	P;
	mat	=	dmatrix(1,M,1,N);
	if(mat == NULL)
	{
		DEF = false;
		errorlog.Add("Matrix::make ","Memory allocation failed");
		return false;
	}
	DEF	=	true;
		for(i=1;i<=M;i++)
			for(j=1;j<=N;j++)
				mat[i][j] = 0.0;
	return DEF;
}

Matrix::~Matrix()
{
	if(DEF)
		free_matrix(mat,1,M,1,N);
}//~Matrix

Matrix& Matrix::operator =(Matrix& B)		//define mat=B
{	
	int i,j;

	if(M!=B.M || N!=B.N || !DEF || !B.DEF)
	{
		errorlog.Add("Matrix::= ","Matrices not compatible for this operation");
		return *this;
	}

	for(i=1;i<=M;i++)
		for(j=1;j<=N;j++)
			mat[i][j] = B(i,j);
		return *this;

}//operator =

Matrix& Matrix::operator +(Matrix& B)		//define A+B
{	
	int i,j;

	if(M!=B.M || N!=B.N || !DEF || !B.DEF)
	{
		errorlog.Add("Matrix::+ ","Matrices not compatible for this operation");
		return *this;
	}
	for(i=1;i<=M;i++)
	{
		for(j=1;j<=N;j++)
			mat[i][j] += B(i,j);
	}//for i
	return *this;
}//operator +

Matrix& Matrix::operator -(Matrix& B)		//define A-B
{	
	int i,j;

	if(M!=B.M || N!=B.N || !DEF)
	{
		errorlog.Add("Matrix::- ","Matrices not compatible for this operation");
		return *this;
	}
	for(i=1;i<=M;i++)
	{
		for(j=1;j<=N;j++)
			mat[i][j] -= B(i,j);
	}//for i
	return *this;
}//operator -

Matrix&	Matrix::operator /(double s)
{
	int i,j;
	for(i=1;i<=M;i++)
	{
		for(j=1;j<=N;j++)
			mat[i][j] /= s;
	}
	return *this;
}

Matrix& Matrix::operator *(double s)
{
	int i,j;
	for(i=1;i<=M;i++)
	{
		for(j=1;j<=N;j++)
			mat[i][j] *= s;
	}
	return *this;
}

Matrix& Matrix::operator *(Matrix& B)		//define A*B
{	
	int i,j,k;
	double temp=0.0;
	Matrix	C(M,N);

	C	=	*this;
	if(N!=B.M || !DEF || !B.DEF)
	{
		errorlog.Add("Matrix::* ","Matrices not compatible for this operation");
		return *this;
	}
	for(i=1;i<=M;i++)
	{
		for(j=1;j<=N;j++)
		{
			temp=0.0;
			for(k=1;k<=N;k++)
				temp+=C(i,k)*B(k,j);
			mat[i][j] = temp;
		}//for j
	}//for i
	return *this;
}//operator*

Vector& Matrix::getrow(int a, Vector& A)
{
	Vector temp;

	if(a<1 || a>M)
	{
		errorlog.Add("Matrix::getrow ","Row index out of bounds");
		return A;
	}
	temp.make(N);
	if(!temp.ok())
		return A;

	for(int i=1;i<=N;i++)
		temp(i,mat[a][i]);
	A = temp;
	return A;
}

Vector& Matrix::getcol(int a, Vector& A)
{
	Vector temp;
	if(a<0||a>N)
	{
		errorlog.Add("Matrix::getcol ","Column index out of bounds");
		return A;
	}
	temp.make(M);
	if(!temp.ok())
		return A;

	for(int i=1;i<=M;i++)
		temp(i,mat[i][a]);
	A = temp;
	return A;
}

Matrix& Matrix::putrow(int a, Vector& A)
{
	if(a<0||a>M)
	{
		errorlog.Add("Matrix::putrow ","Row index out of bounds");
		return *this;
	}
	for(int i=1;i<=N;i++)
		mat[a][i] = A(i);
	return *this;
}

Matrix& Matrix::putcol(int a, Vector& A)
{
	if(a<0||a>N)
	{
		errorlog.Add("Matrix::putcol ","Column index out of bounds");
		return *this;
	}
	for(int i=1;i<=M;i++)
		mat[i][a] = A(i);
	return *this;
}

int	Matrix::rows()			//return the number of rows
{
	if(!DEF)
	{
		errorlog.Add("Matrix::rows ","Matrix not defined");
		return 0;
	}
	return M;
}//rows

int Matrix::cols()			//return the number of columns
{
	if(!DEF)
	{
		errorlog.Add("Matrix::rows ","Matrix not defined");
		return 0;
	}
	return N;

}//cols

bool	Matrix::ok()		{	return	DEF;		}//defined
void Matrix::reset()		
{	
	int i,j;
	for(i=1;i<=M;i++)
		for(j=1;j<=N;j++)
			mat[i][j] = 0.0;
}//reset

double	Matrix::norm2()
{
	int i,j;
	double norm=0.0;
	for(i=1;i<=M;i++)
		for(j=1;j<N;j++)
			norm += mat[i][j]*mat[i][j];

	return sqrt(norm);
}


double	Matrix::operator ()(int i, int j)
{	
	if(i<=0||i>M||j<=0||j>N)
	{
		errorlog.Add("Matrix::(i,j) ","Index out of bounds.");
		return 0.0;
	}
	return mat[i][j];
}

bool	Matrix::operator() (int i, int j, double value)
{	
	if(i<=0||i>M||j<=0||j>N)
	{
		errorlog.Add("Matrix::(i,j,val) ","Index out of bounds.");
		return false;
	}
	mat[i][j] = value;	
	return true; 
}

double	Matrix::fabsrowsum(int i)
{
	double temp=0.0;
	for(int j=1;j<=N;j++)
		temp += fabs(mat[i][j]);
	return temp;
}

double	Matrix::colsum(int j)
{
	double temp=0.0;
	for(int i=1;i<=M;i++)
		temp += fabs(mat[i][j]);
	return temp;
}

double	Matrix::RowInfinityNorm()
{
	double temp1 = 0.0, temp2;
	for(int i=1;i<=M;i++)
	{
		temp2 = this->fabsrowsum(i);
		temp1 = (temp1>temp2)?temp1:temp2;
	}
	return temp1;
}


int		Matrix::ExpSeriesLen()
{
	//Maximum number of terms required to 
	//compute the exponential of the matrix
	//...method of Seem et al (1989), Atkinson (1978)
	int temp;
	temp = int(3*(this->RowInfinityNorm())+6);

	return (temp<100)?temp:100;
}
