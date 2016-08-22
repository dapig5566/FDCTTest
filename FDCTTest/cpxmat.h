#ifndef CPXMAT_H
#define CPXMAT_H
#include "mycomplex.h"
#include <ctime>
class Mat{
public:
    int r,c;
    complex *m;
    Mat(int ri,int ci);
    Mat();
	Mat(const Mat &x);
   ~Mat(); 
      complex& operator ()(int,int);
      void     operator =(Mat );
    friend Mat operator +(Mat , Mat );
    friend Mat operator -(Mat , Mat );
    friend Mat operator *(Mat , Mat );

    friend Mat operator +(Mat , complex );
    friend Mat operator -(Mat , complex );
    friend Mat operator *(Mat , complex );
    friend Mat operator /(Mat , double);
    Mat Multiply(Mat b);
    Mat Cut(int,int,int,int);
    Mat sqrt();
    void Transpose();
	void Reset(int ri, int ci);
};
complex **GetArray(int,int);
void releaseArray(complex**, int);
#endif // CPXMAT_H
