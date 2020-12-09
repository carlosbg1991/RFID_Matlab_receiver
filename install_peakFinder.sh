#!/bin/bash

echo "[INFO] Installing Matlab Peak Finder..."

# download 
(wget -O PeakFinder_Matlab.zip https://terpconnect.umd.edu/~toh/spectrum/PeakFinder.zip --no-check-certificate)

# extract and configure
(unzip PeakFinder_Matlab.zip -d utils/)

# delete the tar file
(rm PeakFinder_Matlab.zip)

echo "[INFO] Succesfully installed Matlab Peak Finder"
echo "[INFO] For more information, visit https://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm"
