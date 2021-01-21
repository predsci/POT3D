#!/bin/bash

rm testsuite/validation/run/* 2>/dev/null
rm bin/pot3d 2>/dev/null

cd src/
make clean 2>/dev/null
rm Makefile 2>/dev/null







