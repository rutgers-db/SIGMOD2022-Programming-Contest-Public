sed -f sed_rule1.txt X1.csv > clean_X1.csv
sed -f sed_rule2.txt X2.csv > temp_X2.csv
sed -f language_rule.txt temp_X2.csv > clean_X2.csv
