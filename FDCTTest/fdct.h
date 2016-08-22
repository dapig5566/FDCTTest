#ifndef FDCT_H
#define FDCT_H
#include "fftw3.h"
#include "cpxmat.h"
#include <iostream>
//#pragma comment(lib,"libfftw3-3.lib")
//#pragma comment(lib,"libfftw3f-3.lib")
//#pragma comment(lib,"libfftw3l-3.lib")
complex**** SepAngles(Mat&, int, int, int&);
Mat** FDCT(Mat& img);
Mat* SepScale(Mat& f, int cor, int J);
Mat* DetailCoeff(Mat& s);
Mat* toClockWise(complex**** R, int S2, int S3, int S4);
Mat Evaluate(complex**f, int t, Mat& shift, int L, Mat W, int n, fftw_complex *in, fftw_complex *out, fftw_plan p);

void fftshift(Mat &m);
void fftshift(complex **&m, int r, int c);
void toMat(complex** R, Mat &f);
void toArray(complex** R, Mat &f);
void LoadMat(Mat& f, fftw_complex*& in);
void OutPut(fftw_complex* out, Mat& f);
void Reshape(complex ****R, int l1, int l2, Mat &t, int L2, int Angles);
void InitFFT(fftw_complex *&in, fftw_complex *&out, int n0, int n1 = 1);
int* GetPoints(int a, int b);

complex**** Squeeze(complex**** R, int L, int LMax, int Angles, int n4);
complex** FFTalongY(Mat s, fftw_complex* in, fftw_complex* out, fftw_plan p);

Mat GetWnd(int l, int sign);
Mat toCol(Mat f);
Mat GetMat(Mat f, int r, int c);
Mat mRound(Mat f);
Mat sum(Mat& Rows, complex **F2, Mat& DeltaP);
Mat PadZero(Mat& t, int n2);
Mat GetSinWnd(int L, double Alpha);
#endif // FDCT_H
