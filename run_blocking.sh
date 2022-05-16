rm output2.csv
rm output1.csv

./blocking clean_X1.csv output1.csv 1000000 0.8 0.57 9 15 & ./blocking clean_X2.csv output2.csv 2000000 0.9 0.57 7 25 &

wait
echo "done blocking"
cat > output.csv <<  EOF
left_instance_id,right_instance_id
EOF
cat output1.csv >> output.csv
cat output2.csv >> output.csv
