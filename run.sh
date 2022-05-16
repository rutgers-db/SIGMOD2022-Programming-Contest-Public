rm output.csv
rm submission.rpz
rm clean_X1.csv
rm clean_X2.csv
rm temp_X2.csv
make clean


reprozip trace --overwrite make
reprozip trace --continue bash sed_clean.sh
reprozip trace --continue bash run_blocking.sh
reprozip pack submission.rpz
