#include <stdio.h>

#define endian endian_

void endian (int*);

void endian (int *endian)
{
    short int one =1;
    *endian = ((int) *(char *) &one);
}
