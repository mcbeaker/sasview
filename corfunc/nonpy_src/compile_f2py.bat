echo off
REM  COMPILE_F2PY

REM  sig.f removed from f2py compilation

REM --fcompiler=gnu95

echo on

gcc -c *.c

python C:\Python27\Scripts\f2py.py -c --fcompiler=gfortran --compiler=mingw32 -lmsvcr90 -m sasview_corfunc addstr.f corfunc_io.f covsrt.f getpath.f mrqmin.f savgol.f strpath.f syminv.f tbksb.f tqli.f tred2.f loqread.f tropus_io.f polint.f sasview_corfunc.f

pause

echo off

REM  Finally tidy up
del *.o
 
pause
