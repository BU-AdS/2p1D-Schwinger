#!/bin/bash

(cd staggered/2D; make clean; rm -rf gauge data logs;)
(cd staggered/3D; make clean; rm -rf gauge data logs;)
(cd wilson/2D; make clean; rm -rf gauge data logs;)
(cd wilson/3D; make clean; rm -rf gauge data logs;)
