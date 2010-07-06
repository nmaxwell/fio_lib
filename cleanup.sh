#!/bin/sh

find . -name *DS_Store* -exec rm {} \;

find . -name *.out -exec rm {} \;
find . -name *.o -exec rm {} \;
find . -name *.a -exec rm {} \;
find . -name *.so -exec rm {} \;
find . -name *.so.1 -exec rm {} \;
find . -name *.so.1.0 -exec rm {} \;

find . -name *.pyc -exec rm {} \;

