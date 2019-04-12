#!/bin/bash

(cd staggered/2D; make clean; rm -rf gauge data logs Makefile main.cpp 2D-Staggered* *~;)
(cd staggered/3D; make clean; rm -rf gauge data logs Makefile main.cpp 2p1D-Staggered *~;)
(cd wilson/2D; make clean; rm -rf gauge data logs Makefile main.cpp 2D-Wilson* *~;)
(cd wilson/3D; make clean; rm -rf gauge data logs Makefile main.cpp 2p1D-Wilson* *~;)
