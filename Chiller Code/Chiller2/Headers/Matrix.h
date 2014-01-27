#ifndef _MATRIX
#define _MATRIX

class	Vector
{
	double *vec;
	int		N;
	bool	OK;
public:
	Vector();
	Vector(int ia);
	~Vector();

	bool	make(int i);
	Vector& operator = (Vector& V);
	Vector& operator + (Vector& V);
	Vector& operator - (Vector& V);
	double	operator() (int i)	{	return vec[i];	}
	bool	operator() (int i, double value)	{	vec[i] = value;	return true; }

	double	min();	//Returns the minimum value in the vector
	double	max();	//Returns the maximum value in the vector
	double	fmin();	//Returns the absolute minimum value in the vector
	double	fmax();	//Returns the absolute maximum value in the vector
	double	sum();	//Returns the sum of all elements in the vector
	double	norm();	//Returns the Euclidean norm of the vector
	void	reset();
	int		size();
	void	scale(double);	//Multiplies all elements by a constant
	bool	ok();
};

//===============================MATRIX CLASS ======================================

class Matrix
{
protected:
	double	**mat;		//The matrix
	int		M,N;		// M-rows, N-columns 
	bool	DEF;		//defined flag	-	true if defined

public:
	//CONSTRUCTORS
	Matrix();					//default constructor - initializes mat to a 3x3 zero matrix
	Matrix(int P);				//alternate constructor - initializes mat to a PxP zero matrix
	Matrix(int P,int Q);		//alternate constructor - initializes mat to a PxQ zero matrix

	bool make(int P, int Q);	//specify size of the matrix after construction

	//DESTRUCTOR
	~Matrix();

	//elementary algebraic operations with matrices
	
	Matrix& operator =(Matrix& B);		//define mat=B
	Matrix& operator +(Matrix& B);		//define A+B
	Matrix& operator -(Matrix& B);		//define A-B
	Matrix& operator /(double s);		//divide each element of matrix by s
	Matrix& operator *(Matrix& B);		//define A*B
	Matrix& operator *(double s);		//multiplies each element of matrix by s
	double	operator() (int i, int j);
	bool	operator() (int i, int j, double value);

	//accessing elements of the matrix and the whole matrix itself
	Vector&	getrow(int i, Vector& A);	//returns the ith row
	Vector& getcol(int i, Vector& A);	//returns the ith column
	Matrix& putrow(int i,Vector& A);	//insert A as the ith row
	Matrix& putcol(int i,Vector& A);	//insert A as the ith column
	int		rows();						//return the number of rows
	int		cols();						//return the number of columns
	bool	ok();

	double	fabsrowsum(int i);			//returns the sum of ith row elements
	double	colsum(int i);				//returns the sum of ith col elements
	double	norm2();					//return the 2-norm of the matrix
	double	RowInfinityNorm();			//returns the row infinity norm of the matrix
	int		ExpSeriesLen();				//returns the number of terms required in the exponential series
										//Seem et al (1989), Atkinson (1978).
	void reset();

};//Class Matrix

//===================================MATRIX CLASS ========================================


#endif