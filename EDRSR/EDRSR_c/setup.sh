#! /bin/sh
swig -python fillvolume.i
gcc -fPIC -c fillvolume.c fillvolume_wrap.c -I /usr/include/python2.7
ld -shared fillvolume.o fillvolume_wrap.o -o _fillvolume.so
