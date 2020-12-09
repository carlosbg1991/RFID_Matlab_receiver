#!/bin/bash

echo "[INFO] Creating required directories"
([ -d misc ] || mkdir misc )
(cd misc ; ([ -d data ] || mkdir data ))

echo "[INFO] Checking on Sample Data set"
(cd misc/data ; ([ -f file_source_test ] || wget https://github.com/nkargas/Gen2-UHF-RFID-Reader/blob/master/gr-rfid/misc/data/file_source_test))

echo "[INFO] Data set is configured correctly"