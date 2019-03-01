#!/bin/bash

(cd staggered/2D; make clean; rm -rf gauge data;)
(cd staggered/3D; make clean; rm -rf gauge data;)
(cd wilson/2D; make clean; rm -rf gauge data;)
(cd wilson/3D; make clean; rm -rf gauge data;)
