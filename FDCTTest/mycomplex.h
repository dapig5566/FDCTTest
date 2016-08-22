#ifndef MYCOMPLEX_H
#define MYCOMPLEX_H
#include <iostream>
class complex{
public:
    double re;       //x
    double im;       //iy
    complex(){re=0;im=0;}
    complex(double a1,double a2=0){re=a1;im=a2;}
    friend complex operator +(complex ,complex);
    friend complex operator -(complex ,complex);
    friend complex operator *(complex ,complex);
    friend complex operator /(complex ,double);
};

#endif // MYCOMPLEX_H
