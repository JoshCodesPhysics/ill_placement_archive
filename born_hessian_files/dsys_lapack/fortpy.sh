#!/bin/bash
python3 disp_solve.py
gfortran lapack_solve.f03 -L/usr/lib -llapack -L/usr/lib -lblas
./a.out
