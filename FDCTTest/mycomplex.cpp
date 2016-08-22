#include "mycomplex.h"

complex operator +(complex a,complex b)
{
    return complex(a.re+b.re,a.im+b.im);
}
complex operator -(complex a,complex b)
{
    return complex(a.re-b.re,a.im-b.im);
}
complex operator *(complex a,complex b)
{
    return complex(a.re*b.re-a.im*b.im , a.re*b.im+a.im*b.re);
}
complex operator /(complex a,double b)
{
    return complex(a.re/b,a.im/b);
}
