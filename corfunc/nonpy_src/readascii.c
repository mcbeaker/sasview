#include <stdio.h>
#include <string.h>

#define readascii readascii_
#define Str_Desc char

void readascii (Str_Desc*,int*,float*,int*);

void readascii (Str_Desc *fname, int *notok, float *data, int* rtrn)
{
    int i,j,iflag,nChannels,nFrames;
    FILE *fp;
    char buff[1000000];
    char string[30];
/* redimensioned tmpdata from 512*512 to match MaxDim in main code, S King, Nov 05 */
    float tmpdata[4096*4096];
    char *buffptr;
    float *dptr;

    i=0;
    do {
        if(fname[i]==' ')
        {
            fname[i]='\0';
            break;
        }
        i++;
    }while(fname[i]);

    if((fp=fopen(fname,"r"))==NULL) {
        *rtrn=0;
        return;
    }

    dptr=tmpdata;
    nChannels=0;
    nFrames=0;
    while(fgets(buff,1000000,fp)!=NULL) {
        buffptr=buff;
        iflag=0;
        nChannels++;
        j=0;
        for(i=0;i<strlen(buff);i++) {
            if(buff[i]==',') {
/* &string has been changed to string in next line as the name string is already a pointer to the 
 first element - removes a compile time warning under Cygwin, S King, Jul 04 */
                strncpy(string,buffptr,i-iflag);
                string[i-iflag]='\0';
                iflag=i+1;
                buffptr=&buff[iflag];
                sscanf(string,"%f",dptr);
                dptr++;
                j++;
            }
        }
        if(nFrames==0)nFrames=j;
        else if(j!=nFrames) {
            *rtrn=0;
            return;
        }
    }

    if(feof(fp)==0||ferror(fp)!=0) {
        *rtrn=0;
        return;
    }

    notok[0]=nChannels;
    notok[1]=nFrames;
    notok[2]=1;
    for(i=3;i<10;i++) {
        notok[i]=0;
    }

    dptr=data;
    for(j=0;j<nFrames;j++) {
        for(i=0;i<nChannels;i++) {
            *dptr=tmpdata[(i*nFrames)+j];
            dptr++;
        }
    }

    fclose(fp);
    *rtrn=1;
    return;
}
