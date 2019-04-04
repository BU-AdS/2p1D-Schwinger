#!/bin/bash

(cd staggered/2D; make clean; rm -rf gauge data logs Makefile main.cpp;)
(cd staggered/3D; make clean; rm -rf gauge data logs Makefile main.cpp;)
(cd wilson/2D; make clean; rm -rf gauge data logs Makefile main.cpp 2D-Wilson*;)
(cd wilson/3D; make clean; rm -rf gauge data logs Makefile main.cpp 2D-Wilson*;)
