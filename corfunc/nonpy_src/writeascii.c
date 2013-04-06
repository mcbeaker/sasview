#include <stdio.h>

#define writeascii writeascii_
#define Str_Desc char

void writeascii (Str_Desc*,int*,float*,int*,int*);

void writeascii (Str_Desc *fname, int *notok, float *data, int* rtrn, int *flag)
{
    int nChannels,nFrames,i,j;
    FILE *fp;
    float *dptr;

    nChannels=notok[0];
    nFrames=notok[1];

    i=0;
    do {
        if(fname[i]==' ')
        {
            fname[i]='\0';
            break;
        }
        i++;
    }while(fname[i]);

    if(flag) {
        if((fp=fopen(fname,"w"))==NULL) {
            *rtrn=0;
            return;
        }
    } else {
        if((fp=fopen(fname,"a"))==NULL) {
            *rtrn=0;
            return;
        }
    }

    for(i=0;i<nChannels;i++) {
        dptr=data+i;
        for(j=0;j<nFrames;j++) {
            fprintf(fp,"  %7.5E,",*dptr);
            dptr+=nChannels;
        }
        fprintf(fp,"\n");
    }

    if(ferror(fp)!=0) {
        *rtrn=0;
        return;
    }

    fclose(fp);
    *rtrn=1;
    return;
}
