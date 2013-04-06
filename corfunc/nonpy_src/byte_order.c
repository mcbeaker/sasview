#include <stdio.h>

#define swap swap_
#define endian endian_

void swap (void*,int*,int*,int*);
void swabytes (void**, int, int);
void endian (int*);

void endian (int* nval)
{
    short int one =1;
    *nval = ((int) *(char *) &one);
}

void swap (void *vp, int* fendian, int* n, int* dsize)
{
    int swab,i,nval;
    endian(&nval);

    swab=*fendian^nval;

    if(swab)
      swabytes(&vp,*n,*dsize);
}

void swabytes (void **vp, int n, int dsize)
{
    unsigned char *cp = (unsigned char *)*vp;
    int t,i,j;

    while (n-- > 0)
    {
       j=0;
       for(i=dsize-1;i>j;i--)
       {
         t = cp[i]; cp[i] = cp[j]; cp[j] = t;
         j++;
       }
       cp+=dsize;
    }
}
