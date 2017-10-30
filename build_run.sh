#!/bin/bash
cd /mnt/d/Studies-MOOCs/2017-02_Self-Driving-Car_NanoDegree_Udacity/Self_Driving_Car_Projects/9_Path_Planning_Project/CarND-Path-Planning-Project
rm -r build
mkdir build && cd build
cmake ..
make
./path_planning
