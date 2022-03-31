#!/bin/bash

rm testsuite/validation/run/* 2>/dev/null
touch testsuite/validation/run/empty
rm bin/pot3d 2>/dev/null

cd src/
make clean 1>/dev/null 2>/dev/null
rm build.log 2>/dev/null
rm build.err 2>/dev/null
rm Makefile 2>/dev/null







