rm output.csv
rm submission.rpz
rm clean_X1.csv
rm clean_X2.csv
rm temp_X2.csv
make clean


make
bash sed_clean.sh
bash run_blocking.sh

