#!/bin/bash

rm testsuite/validation/run/* 2>/dev/null
rm testsuite/isc2023/run/* 2>/dev/null
rm testsuite/small/run/* 2>/dev/null
rm testsuite/medium/run/* 2>/dev/null
rm testsuite/large/run/* 2>/dev/null
rm bin/pot3d 2>/dev/null

cd src/
make clean 1>/dev/null 2>/dev/null
rm build.log 2>/dev/null
rm build.err 2>/dev/null
rm Makefile 2>/dev/null







