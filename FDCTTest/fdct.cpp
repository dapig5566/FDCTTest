#include "fdct.h"
#include <cmath>
#define COARSEST 4
#define pi 3.1415926

using std::cout;
using std::endl;
Mat** FDCT(Mat& img)
{
    fftw_plan p=NULL;
	fftw_complex *in = NULL, *out = NULL;
    int J=log2(img.r*1.0);
    Mat** CM=new Mat*[J-COARSEST+1];
	time_t s, e;
	s = clock();
    Mat* Scales=SepScale(img,COARSEST,J);
	e = clock();
	cout << (e-s)*1.0 / CLOCKS_PER_SEC << endl;
    for(int i=0;i<J-COARSEST+1;i++)
    {
        if(i==0||i==J-COARSEST)
        {

            InitFFT(in,out,Scales[i+1].r,Scales[i+1].c);
            p=fftw_plan_dft_2d(Scales[i+1].r,Scales[i+1].c,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
            s = clock();
            fftshift(Scales[i+1]);
            LoadMat(Scales[i+1],in);
            fftw_execute(p);
            CM[i]=new Mat(Scales[i+1].r,Scales[i+1].c);
            OutPut(out, CM[i][0]);
            fftshift(CM[i][0]);
            CM[i][0]=CM[i][0]/::sqrt(Scales[i+1].r*Scales[i+1].c*1.0);
            e = clock();
            cout << (e - s)*1.0 / CLOCKS_PER_SEC << endl;
			fftw_free(in);
            fftw_free(out);
			fftw_destroy_plan(p);
        }else{
            s = clock();
          CM[i]=DetailCoeff(Scales[i+1]);
            e = clock();
            cout << (e - s)*1.0 / CLOCKS_PER_SEC << endl;
        }
    }
	delete[] Scales;
    return CM;
}
Mat* SepScale(Mat& f,int cor,int J)
{
    int n=f.r;
    fftw_plan p;
    fftw_complex *in=NULL,*out=NULL;
    Mat F(f.r,f.c),
        *S=new Mat[J-cor+2];


/*
       for(int i=1;i<=10;i++)
        cout<<f(1,i).re<<"+("<<f(1,i).im<<"i)\n";

*/  InitFFT(in,out,f.r,F.c);
    p=fftw_plan_dft_2d(f.r,f.c,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
    fftshift(f);
    LoadMat(f,in);
    fftw_execute(p);
    OutPut(out,F);
    fftshift(F);
    F=F/::sqrt(f.r*f.c*1.0);
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    for(int i=cor;i<=J;i++)
    {
        Mat Whi1,Whi2,WLo1,WLo2;
		int p = pow(2, i);
		int l=p<n/2?p:n/2;
        if(i==cor)
        {
            Whi1=GetWnd(i,0);
            Whi2=Whi1;
            Whi2.Transpose();
            Whi2=Whi2*Whi1;
        }else{
            WLo1=GetWnd(i-1,1);
            WLo2=WLo1;
            WLo2.Transpose();
            WLo2=WLo2*WLo1;
			WLo2 = WLo2.Multiply(WLo2);
            if(i==J)
            {
				Whi1.Reset(1, 2 * l);
                for(int i=1;i<=2*l;i++)
                    Whi1(1,i)=complex(1);
            }else{
				Whi1=GetWnd(i,0);
            }
			Whi2 = Whi1;
			Whi2.Transpose();
			Whi2 = Whi2*Whi1;
			Whi2 = Whi2.Multiply(Whi2);
			Whi2 = Whi2 - WLo2;
			Whi2 = Whi2.sqrt();
        }
		S[i - cor + 1] = Whi2.Multiply(F.Cut(n / 2 - l + 1, n / 2 + l, n / 2 - l + 1, n / 2 + l));
    }
    return S;
}
void fftshift(Mat &m)
{
    Mat tmp(m.r,m.c);
	int cnt = 1,
		*x=new int[m.r+1],*y=new int[m.c+1];
    for(int i=1+m.r/2;i<=m.r;i++)
        x[i]=cnt++;
    for(int i=1;i<=m.r/2;i++)
        x[i]=cnt++;
    cnt=1;
    for(int i=1+m.c/2;i<=m.c;i++)
        y[i]=cnt++;
    for(int i=1;i<=m.c/2;i++)
        y[i]=cnt++;

    for(int i=1;i<=m.r;i++)
        for(int j=1;j<=m.c;j++)
            tmp(i,j)=m(x[i],y[j]);
    m=tmp;
    delete[] x;
    delete[] y;
}
void fftshift(complex **&m,int r,int c)
{
    complex **tmp=GetArray(r,c);
    int cnt = 1,
        *x=new int[r+1],*y=new int[c+1];
    for(int i=1+r/2;i<=r;i++)
        x[i]=cnt++;
    for(int i=1;i<=r/2;i++)
        x[i]=cnt++;
    cnt=1;
    for(int i=1+c/2;i<=c;i++)
        y[i]=cnt++;
    for(int i=1;i<=c/2;i++)
        y[i]=cnt++;

    for(int i=1;i<=r;i++)
        for(int j=1;j<=c;j++)
            tmp[i][j]=m[x[i]][y[j]];
    releaseArray(m,r);
    m=tmp;
    delete[] x;
    delete[] y;
}
Mat* DetailCoeff(Mat& s)
{
    int Angles=pow(2,floor((log2(s.r/2)+1)/2)),
        L=s.r/Angles,
        S1=4,S2=Angles,S3=s.r/4,S4=0;
    int cnt=1,
        *x,*y;
    complex ****R;
	fftw_complex *in = NULL, *out = NULL;
    fftw_plan p=NULL;

    R=SepAngles(s,Angles,L,S4);
	x = new int[S3+1];
	y = new int[S4+1];
    for(int i=2;i<=S3;i++)
        x[i]=cnt++;
    x[1]=S3;
    cnt=1;
    for(int i=2;i<=S4;i++)
        y[i]=cnt++;
    y[1]=S4;

    for(int k=1;k<=S2;k++)
    {
        Mat tmp1(S3,S4),tmp2(S3,S4);
        for(int i=1;i<=S3;i++)
            for(int j=1;j<=S4;j++)
            {
                tmp1(i,j)=R[2][k][x[i]][y[j]];
                tmp2(i,j)=R[4][k][x[i]][y[j]];
            }
        for(int i=1;i<=S3;i++)
            for(int j=1;j<=S4;j++)
            {
				R[2][k][i][j] = tmp1(i, j);
				R[4][k][i][j] = tmp2(i, j);
            }
    }
    delete[] x;
    delete[] y;
    InitFFT(in,out,S3,S4);
	p = fftw_plan_dft_2d(S3, S4, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    for(int i=1;i<=S1;i++)
        for(int j=1;j<=S2;j++)
        {
            Mat tmp(S3,S4);
            toMat(R[i][j],tmp);
            fftshift(tmp);
            LoadMat(tmp,in);
            fftw_execute(p);
            OutPut(out,tmp);
			fftshift(tmp);
            tmp=tmp/::sqrt(S3*S4*1.0);
            toArray(R[i][j],tmp);
        }
    fftw_free(in);
    fftw_free(out);
    fftw_destroy_plan(p);
    return toClockWise(R,S2,S3,S4);

}
complex**** SepAngles(Mat& S,int Angles,int L,int& S4)
{
    int n=S.r,LMax,
        *Index;
    Mat       M(1,Angles),
              Y,Slope;
    complex **F1,**F2,****R;
	fftw_complex *in, *out, *inT, *outT;
	fftw_plan p,pT;
    double AlphaMax;
    Index=GetPoints(n/8,n/4);

    R=new complex***[5];
    for(int i=1;i<=4;i++)
        R[i]=new complex**[Angles+1];
    for(int i=1;i<=4;i++)
        for(int j=1;j<=Angles;j++)
            R[i][j]=GetArray(S.r/4,2*L);
    
	
	for(int i=1;i<=n/8;i++)
    {
        int tmp;
        tmp = -Index[i];
        Index[i] = -Index[n/4-i+1];
        Index[n/4-i+1] = tmp;
    }
	
    for(int i=1;i<=M.c;i++)
        M(1,i)=complex(i-1);
	
	Y=(M-complex(Angles*1.0/2))*complex(L)+complex(L*1.0/2);
	
    Slope=Y;
    for(int i=1;i<=Slope.c;i++)
		Slope(1, i).re *=- 1;

    Slope=Slope/(n*1.0/2);
	InitFFT(in, out, S.r);
	p = fftw_plan_dft_1d(S.r, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    F1=FFTalongY(S,in,out,p);
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
	for(int i=1;i<=S.r;i++)
        for(int j=1;j<=S.c;j++)
            F1[i][j]=F1[i][j]/::sqrt(S.r*1.0);
    S.Transpose();
	InitFFT(in, out, S.r);
	p = fftw_plan_dft_1d(S.r, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    F2=FFTalongY(S,in,out,p);
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
    for(int i=1;i<=S.r;i++)
        for(int j=1;j<=S.c;j++)
            F2[i][j]=F2[i][j]/::sqrt(S.r*1.0);

	if (Angles <= 16)
	{
		InitFFT(in, out, S.c);
		InitFFT(inT, outT, S.r);
		p = fftw_plan_dft_1d(S.c, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		pT = fftw_plan_dft_1d(S.r, inT, outT, FFTW_FORWARD, FFTW_ESTIMATE);
	}else{
		InitFFT(in, out, S.c * 16);
		InitFFT(inT, outT, S.r * 16);
		p = fftw_plan_dft_1d(S.c * 16, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
		pT = fftw_plan_dft_1d(S.r * 16, inT, outT, FFTW_FORWARD, FFTW_ESTIMATE);
	}
    for(int i=1;i<=n/4;i++)
    {
		int		t = Index[i] + n / 2 + 2;
        double  Alpha;
        Mat Shift,
            W,tmp;
        Shift=Y+Slope*complex(t-1);
        Alpha=-(Index[i]+1)*1.0/(n/2);

        W=GetSinWnd(L,Alpha);
		tmp = Evaluate(F1, t, Shift, L, W, S.c, in, out, p);
        Reshape(R,1,i,tmp,2*L,Angles);

		tmp = Evaluate(F2, t, Shift, L, W, S.r, inT, outT, pT);
        Reshape(R,3,i,tmp,2*L,Angles);


        t=n/2-Index[i];
        for(int i=1;i<=Shift.c;i++)
            Shift(1,i).re*=-1;
        Shift=Shift+complex(1);

		tmp = Evaluate(F1, t, Shift, L, W, S.c, in, out, p);
        Reshape(R,2,n/4+1-i,tmp,2*L,Angles);

		tmp = Evaluate(F2, t, Shift, L, W, S.r, inT, outT, pT);
        Reshape(R,4,n/4+1-i,tmp,2*L,Angles);

    }
	fftw_destroy_plan(p);
	fftw_destroy_plan(pT);
	fftw_free(in);
	fftw_free(out);
	fftw_free(inT);
	fftw_free(outT);
    delete[] Index;
    releaseArray(F1,S.r);
    releaseArray(F2,S.r);
    Index=GetPoints(n/8,n/4);
    AlphaMax=(Index[n/4]-1)*1.0/(n/2);
    LMax=ceil(AlphaMax*L*1.0);
    S4=2*LMax;
	delete[] Index;
    return Squeeze(R,L,LMax,Angles,n/4);
}
Mat* toClockWise(complex**** R,int S2,int S3,int S4)
{
    Mat *D=new Mat[4*S2];
    int cnt=0;
    for (int i = 0; i < S2; i++)
		D[i].Reset(S3, S4);
	for (int i = S2; i < 2 * S2; i++)
		D[i].Reset(S4, S3);
	for (int i = 2*S2; i < 3 * S2; i++)
		D[i].Reset(S3, S4);
	for (int i = 3*S2; i < 4 * S2; i++)
		D[i].Reset(S4, S3);
    for(int i=1;i<=S2;i++)
    {
        for(int j=1;j<=S3;j++)
            for(int k=1;k<=S4;k++)
                D[cnt](j,k) = R[3][i][j][k];
		cnt++;
    }
    for(int i=S2;i>=1;i--)
    {
        for(int j=1;j<=S3;j++)
            for(int k=1;k<=S4;k++)
                D[cnt](k,j) = R[2][i][j][k];
		
		cnt++;
    }
	
    for(int i=1;i<=S2;i++)
    {
        for(int j=1;j<=S3;j++)
            for(int k=1;k<=S4;k++)
				D[cnt](j, k) = R[4][i][j][k];
		cnt++;
    }
    for(int i=S2;i>=1;i--)
    {
        for(int j=1;j<=S3;j++)
            for(int k=1;k<=S4;k++)
				D[cnt](k, j) = R[1][i][j][k];
		
		cnt++;
    }
	
	for (int i = 1; i <= 4; i++)
		for (int j = 1; j <= S2; j++)
			releaseArray(R[i][j], S3);
	for (int i = 1; i <= 4; i++)
		delete[] R[i];
	delete[] R;
    return D;
}
Mat Evaluate(complex**f, int t, Mat& shift, int L, Mat W, int n, fftw_complex *in, fftw_complex *out, fftw_plan p)
{
    int Angles=shift.c;
    Mat F(n,1),T(n,1),
        E;
	complex **F2;
    for(int i=1;i<=n;i++)
		T(i, 1) = complex(-n*1.0 / 2 + i - 1);
    for(int i=1;i<=n;i++)
        F(i,1)=f[i][t];

    if(Angles<=16)
    {
        Mat Tmp(n,Angles);
        int *Col=new int[2*L+1];
        F=GetMat(F,1,Angles);
        for(int i=1;i<=2*L;i++)
            Col[i]=n/2-L+i;
		Tmp = T*shift / (n*1.0)*complex(-2 * pi);
        for(int i=1;i<=Tmp.r;i++)
        {
			for (int j = 1; j <= Tmp.c; j++)
			{
				Tmp(i, j).im = sin(Tmp(i, j).re);
				Tmp(i, j).re = cos(Tmp(i, j).re);
			}
        }
        F=F.Multiply(Tmp);
        F2=FFTalongY(F,in,out,p);
        for(int i=1;i<=F.r;i++)
            for(int j=1;j<=F.c;j++)
                F2[i][j]=F2[i][j]/::sqrt(n*1.0);
        W.Transpose();
		W = GetMat(W, 1, Angles);
        for(int i=1;i<=W.r;i++)
            for(int j=1;j<=W.c;j++)
                W(i,j)=W(i,j)*F2[Col[i]][j];
		delete[] Col;
		releaseArray(F2, F.r);
        return toCol(W);
    }else{
        Mat K(1,2*L),Rows(2*L*Angles,1),
            t1(2*L,Angles),t2(Angles,2*L),
            Near(2*L*Angles,1),
            Delta(2*L*Angles,1),DeltaP(2*L*Angles,6);

        for(int i=1;i<=2*L;i++)
            K(1,i)=-L+i-1;
        t2=GetMat(K,Angles,1);
        t2.Transpose();
        t2=toCol(t2+GetMat(shift,2*L,1))*complex(2*pi)/n;
        t1=GetMat(F,1,6);
        for(int i=2;i<=6;i++)
            for(int j=1;j<=n;j++)
                t1(j,i)=t1(j,i-1)*complex(0,-T(j,1).re);
        F2=FFTalongY(PadZero(t1,n/2),in,out,p);								//16*n
        Near=mRound(t2*complex(8*n/pi));
        Rows=Near+complex(8*n+1);
        Delta=t2-Near/(8.0*n);
        for(int i=1;i<=6;i++)
            for(int j=1;j<=2*L*Angles;j++)
                DeltaP(j,i)=DeltaP(j,i-1)*Delta(j,1)/(i-1);
        E=sum(Rows,F2,DeltaP);
		releaseArray(F2,15*n+t1.r);
        E=E.Multiply(GetMat(toCol(W),Angles,1));
		return E;
    }
	return Mat();
}

void LoadMat(Mat& f, fftw_complex *&in)
{
    for(int i=0;i<f.r;i++)
        for(int j=0;j<f.c;j++)
        {
            in[i*f.c+j][0]=f(i+1,j+1).re;
            in[i*f.c+j][1]=f(i+1,j+1).im;
        }
}

void OutPut(fftw_complex* out,Mat& f)
{
    for(int i=0;i<f.r;i++)
        for(int j=0;j<f.c;j++)
            f(i+1,j+1)=complex(out[i*f.c+j][0],out[i*f.c+j][1]);
}

void InitFFT(fftw_complex *&in, fftw_complex *&out, int n0, int n1)
{
    in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n0*n1);
    out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n0*n1);
}

Mat GetWnd(int l, int sign)
{
    int dyadic=pow(2,(l-1)),epsp=dyadic/3,
        cnt=2;
    Mat tmp1(1,(2+2*sign)*dyadic),
        tmp2(1,(4+4*sign)*dyadic);
    for(int i=1;i<=dyadic-epsp;i++)
        tmp1(1,i)=complex(1);
    for(int i=dyadic-epsp+1;i<=dyadic+epsp+1;i++)
    {
        double x=(3*1.0/2)*((i-1)*1.0/dyadic)-1;
        if(x<=0)
            x=0;
		else{
			if (x >= 1)
				x = 1;
			else
				x = x*x*x*x*(35 - 84 * x + 70 * x*x - 20 * x*x*x);
		}
        tmp1(1,i)=complex(cos(pi/2*x));
    }
    for(int i=dyadic+epsp+2;i<=2*dyadic;i++)
        tmp1(1,i)=complex(0);

    tmp2(1,1)=complex(0);
    for(int i=tmp1.c;i>=2;i--)
        tmp2(1,cnt++)=tmp1(1,i);
    for(int i=1;i<=tmp1.c;i++)
        tmp2(1,cnt++)=tmp1(1,i);
    return tmp2;
}
int *GetPoints(int a,int b)
{
    int *p=new int[b+1],
        cnt=1,eps=a/3,epsp=a-eps;
    for(int i=epsp+1;i<=b+epsp;i++)
        p[cnt++]=i;
    return p;
}
complex** FFTalongY(Mat s, fftw_complex* in, fftw_complex* out, fftw_plan p)
{
    complex** c=GetArray(s.r,s.c);
    fftshift(s);
    for(int i=1;i<=s.c;i++)
    {
        for(int j=1;j<=s.r;j++)
        {
            in[j-1][0]=s(j,i).re;
            in[j-1][1]=s(j,i).im;
        }
        fftw_execute(p);
        for(int j=1;j<=s.r;j++)
        {
            c[j][i].re=out[j-1][0];
            c[j][i].im=out[j-1][1];
        }
    }
    fftshift(c,s.r,s.c);
    return c;
}
void Reshape(complex**** R,int l1,int l2,Mat &t,int L2,int Angles)
{
    Mat tmp(L2,Angles);
    int x=1,y=1;
    for(int i=1;i<=Angles;i++)
        for(int j=1;j<=L2;j++)
        {
            tmp(j,i)=t(x++,y);
            if(x>t.r)
            {
                x=1;
                y++;
            }
        }
    tmp.Transpose();
    for(int i=1;i<=Angles;i++)
        for(int j=1;j<=L2;j++)
            R[l1][i][l2][j]=tmp(i,j);
}
complex**** Squeeze(complex**** R, int L, int LMax, int Angles, int n4)
{
	int *mid = new int[2 * LMax + 1], cnt = 1;
    complex ****R2;

    R2=new complex***[5];
    for(int i=1;i<=4;i++)
        R2[i]=new complex**[Angles+1];
    for(int i=1;i<=4;i++)
        for(int j=1;j<=Angles;j++)
            R2[i][j]=GetArray(n4,2*LMax);


    for(int i=L-LMax+1;i<=L+LMax;i++)
        mid[cnt++]=i;
    for(int i=1;i<=2*LMax;i++)
        for(int j=1;j<=4;j++)
            for(int k=1;k<=Angles;k++)
                for(int l=1;l<=n4;l++)
                    R2[j][k][l][i]=R[j][k][l][mid[i]];
    delete[] mid;
    for(int i=1;i<=4;i++)
        for(int j=1;j<=Angles;j++)
            releaseArray(R[i][j],n4);
    for(int i=1;i<=4;i++)
        delete[] R[i];
    delete[] R;
    return R2;
}
Mat toCol(Mat f)
{
    Mat tmp(f.r*f.c,1);
    int x=1,y=1;
    for(int i=1;i<=tmp.r;i++)
    {
        tmp(i,1)=f(x++,y);
        if(x>f.r)
        {
            x=1;
            y++;
        }
    }
    return tmp;
}
Mat GetMat(Mat f, int r, int c)
{
    Mat tmp(r*f.r,c*f.c);
    for(int i=1;i<=r*f.r;i++)
        for(int j=1;j<=c;j++)
            memcpy(&tmp.m[i*(c*f.c+1)+(j-1)*f.c+1],&f.m[((i-1)%f.r+1)*(f.c+1)+1],sizeof(complex)*f.c);
	return tmp;
}
Mat mRound(Mat f)
{
	for (int i = 1; i <= f.r;i++)
		for (int j = 1; j <= f.c; j++)
			f(i, j) = complex(floor(f(i, j).re + 0.5));
		return f;
}
Mat  sum(Mat& Rows, complex **F2, Mat& DeltaP)
{
	Mat tmp(Rows.r, DeltaP.c),
		tmp2(Rows.r,1);
	for (int i = 1; i <= tmp.r;i++)
		for (int j = 1; j <= tmp.c; j++)
			tmp(i, j) = F2[(int)Rows(i, 1).re][j] * DeltaP(i, j);
	for (int i = 1; i <= tmp.r; i++)
	{
		complex s(0);
		for (int j = 1; j <= tmp.c; j++)
			s = s + tmp(i, j);
		tmp2(i, 1) = s;
	}
	return tmp2;
}
Mat PadZero(Mat& t, int n2)
{
	Mat tmp(15 * 2 * n2 + t.r, t.c);
	for (int i = 1; i <= 15 * n2;i++)
		for (int j = 1; j <= tmp.c; j++)
			tmp(i, j) = complex(0);
	for (int i = 15 * n2+1; i <= 15 * n2+t.r; i++)
		for (int j = 1; j <= tmp.c; j++)
			tmp(i, j) = t(i - 15 * n2, j);
	for (int i = 15 * n2 + t.r + 1; i <= tmp.r; i++)
		for (int j = 1; j <= tmp.c; j++)
			tmp(i, j) = complex(0);
	return tmp;
}
Mat GetSinWnd(int L, double Alpha)
{
	Mat tmp(1, 2 * L);
	int cnt = -L;
	double *idx = new double[2 * L+1];
	for (int i = 1; i <= 2 * L; i++)
	{
        idx[i] = (cnt + 0.5) / (L*1.0*Alpha);
		cnt++;
	}
	for (int i = 1; i <= 2 * L; i++)
	{
		double  val = 0;
		if (idx[i] <= -1 || idx[i] >= 1)
			val = 0;
		else{
			if (idx[i] == 0)
				val = 1;
			else
				val = idx[i] > 0 ? sin(pi / 4 * (1 + sin(pi*(0.5 - idx[i])))) : sin(pi / 4 * (1 + sin(pi*(0.5 + idx[i]))));
		}
		tmp(1, i) = complex(val);
	}
	return tmp;
}

void toMat(complex **R, Mat &f)
{
	for (int i = 1; i <= f.r;i++)
		for (int j = 1; j <= f.c; j++)
			f(i, j) = R[i][j];
		
}
void toArray(complex **R, Mat &f)
{
	for (int i = 1; i <= f.r; i++)
		for (int j = 1; j <= f.c; j++)
			R[i][j] = f(i, j);
}
