#include "fdct.h"
//#include <ctime>
#include <cstdio>

using namespace std;

int main()
{
	Mat test(512, 512);
	int m = 0;
	int& q = m;
	Mat **C;
	time_t s, e;
	//srand(clock());
	freopen("test.txt", "r",stdin);
	for (int i = 1; i <= 512;i++)
	for (int j = 1; j <= 512; j++)
	{
		int tmp = 0;
		std::cin >> tmp;
		test(i, j) = complex(tmp);
	}
	freopen("CON", "r", stdin);
//	s = clock();
	C=FDCT(test);
//	e = clock();
//	cout << (e - s)*1.0 / CLOCKS_PER_SEC;
/*	Mat t1(4, 4);
	Mat t2(t1),
		test(4,4);
	int cnt = 0;
	for (int i = 1; i <= 4; i++)
		for (int j = 1; j <= 4; j++)
			test(i,j) = complex(cnt++);
	//test = GetMat(t1, 3, 3);
	for (int i = 1; i <= test.r; i++)
	{
		for (int j = 1; j <= test.c; j++)
			cout << test(i, j).re << " ";
		cout << endl;
	}
	cout << endl;
	fftshift(test);
	fftshift(test);
	for (int i = 1; i <= test.r; i++)
	{
		for (int j = 1; j <= test.c; j++)
			cout << test(i, j).re << " ";
		cout << endl;
	}*/
/*	for (int i = 1; i <= 6; i++)
		delete[] C[i];
	delete[] C;*/
	system("pause");
	return 0;
}