#include <stdio.h>
#include <string.h>
#include <ctype.h>

/* Added by SMK, 05/04/13 */
#include <windows.h>


#define getfiletype getfiletype_
#define Str_Desc char

void getfiletype (Str_Desc*,Str_Desc*);

void getfiletype (Str_Desc *filetype, Str_Desc *fname)
{
    int i;
    char *filename;
    char *fname2;

    i=0;
    do {
        if(fname[i]==' ')
        {
            fname[i]='\0';
            break;
        }
        i++;
    }while(fname[i]);

    fname2 = (char*)malloc((strlen(fname)+1)*sizeof(char));

    for(i=0;i<strlen(fname);i++) {
        fname2[i]=(char)(tolower((int)fname[i]));
    }
    fname2[i]='\0';

    for(i=strlen(fname2);i>=0;i--) {
        if(fname2[i]=='.') {
            filename=&fname2[i];
            break;
        }
    }
    if(strcmp(filename,".txt")==0) {
        strcpy(filetype,"ascii     ");}
    else if(strcmp(filename,".dat")==0) {
        strcpy(filetype,"ascii     ");
    } else {
        strcpy(filetype,"nonascii  ");
    }
}
