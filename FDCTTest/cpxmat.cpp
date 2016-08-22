#include "cpxmat.h"
#include <cmath>
Mat::Mat()
{
    r=c=0;
    m=NULL;
}
Mat::Mat(int ri, int ci)
{
    r=ri;
    c=ci;
	m = new complex[(r + 1)*(c + 1)];
}
Mat::~Mat()
{
	delete[] m;
}
complex& Mat::operator ()(int i,int j)
{
    return m[i*(c+1)+j];
}
Mat operator +(Mat a,Mat b)
{
    Mat tmp(a.r,a.c);
    if(a.r==b.r&&a.c==b.c)
    {
		for (int i = a.c + 1+1; i < (a.r + 1)*(a.c + 1); i++)
			tmp.m[i]=a.m[i]+b.m[i];
    }
    return tmp;
}

Mat operator -(Mat a,Mat b)
{
    Mat tmp(a.r,a.c);
    if(a.r==b.r&&a.c==b.c)
    {
		for (int i = a.c + 1 + 1; i < (a.r + 1)*(a.c + 1); i++)
			tmp.m[i] = a.m[i] - b.m[i];
    }
    return tmp;
}

Mat operator *(Mat a,Mat b)
{
    Mat tmp(a.r,b.c);
    if(a.c==b.r)
    {
        for(int i=1;i<=a.r;i++)
            for(int j=1;j<=b.c;j++)
            {
                complex tmp2;
                for(int k=1;k<=a.c;k++)
                    tmp2=tmp2+(a(i,k)*b(k,j));
                tmp(i,j)=tmp2;
            }
    }
    return tmp;
}

Mat operator +(Mat a,complex b)
{
    Mat tmp(a.r,a.c);
	for (int i = a.c + 1 + 1; i < (a.r + 1)*(a.c + 1); i++)
		tmp.m[i] = a.m[i] + b;
    return tmp;
}

Mat operator -(Mat a,complex b)
{
    Mat tmp(a.r,a.c);
	for (int i = a.c + 1 + 1; i < (a.r + 1)*(a.c + 1); i++)
		tmp.m[i] = a.m[i] - b;
    return tmp;
}

Mat operator *(Mat a,complex b)
{
    Mat tmp(a.r,a.c);
	for (int i = a.c + 1 + 1; i < (a.r + 1)*(a.c + 1); i++)
		tmp.m[i] = a.m[i] * b;
    return tmp;
}

Mat operator /(Mat a,double b)
{
    Mat tmp(a.r,a.c);
	for (int i = a.c + 1 + 1; i < (a.r + 1)*(a.c + 1); i++)
		tmp.m[i] = a.m[i] / b;
    return tmp;
}
Mat Mat::Multiply(Mat b)
{
    Mat tmp(r,c);
	if (r == b.r&&c == b.c)
	{
		for (int i = c + 1 + 1; i < (r + 1)*(c + 1); i++)
			tmp.m[i] = m[i] * b.m[i];
	}
    return tmp;
}

Mat Mat::Cut(int a,int b,int c,int d)
{
    Mat tmp(b-a+1,d-c+1);
    int cntr=1,cntc=1;
    for(int i=a;i<=b;i++)
    {
        cntc=1;
        for(int j=c;j<=d;j++)
        {
            tmp(cntr,cntc)=this->operator()(i,j);
            cntc++;
        }
        cntr++;
    }
    return tmp;
}

Mat Mat::sqrt()
{
    Mat tmp(r,c);
	for (int i = c + 1 + 1; i < (r + 1)*(c + 1); i++)
		tmp.m[i] = ::sqrt(m[i].re);
    return tmp;
}
void Mat::Transpose()
{
    complex *t=new complex[(r+1)*(c+1)];
    int tmp;
    for(int i=1;i<=c;i++)
        for(int j=1;j<=r;j++)
            t[i*(r+1)+j]=this->operator()(j,i);
	delete[] m;
	tmp=r;
    r=c;
    c=tmp;
    m=t;
}
void Mat::operator =(Mat b)
{
	r = b.r;
	c = b.c;
	if (m != NULL)
		delete[] m;
	if (b.m == NULL)
		m = NULL;
	else{
		m = new complex[(r + 1)*(c + 1)];
		memcpy(m, b.m, sizeof(complex)*(c + 1)*(r + 1));
	}
}

complex **GetArray(int r,int c)
{
    complex **t=new complex*[r+1];
    for(int i=1;i<=r;i++)
        t[i]=new complex[c+1];
    return t;
}
void releaseArray(complex** R,int r)
{
    for(int i=1;i<=r;i++)
        delete[] R[i];
    delete[] R;
}

void Mat::Reset(int ri, int ci)
{
	r = ri;
	c = ci;
	if (m != NULL)
		delete[] m;
	m = new complex[(r+1)*(c+1)];
}
Mat::Mat(const Mat &x)
{
	m = new complex[(x.r + 1)*(x.c + 1)];
	r = x.r;
	c = x.c;
	memcpy(m, x.m, sizeof(complex)*(r + 1)*(c + 1));
}