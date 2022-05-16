# ACM SIGMOD Programming Contest 2022

## Team Info 

Team Name: April

Member: Chaoji Zuo, Zhizhi Wang

Advisor: [Dong Deng](https://people.cs.rutgers.edu/~dd903/)

Affiliation: 

Data Curation Lab in Database Group,
Computer Science Department at Rutgers University.

## Description

Our approach leverages similarity grouping, similarity join and overlap join to efficiently find a pool of candidate record pairs. After that, we rank all the candidate record pairs by their weighted similarities and output a fixed number of record pairs as the blocking result. We used efficient Jaccard similarity join ([Deng et al. VLDB 2016](https://people.cs.rutgers.edu/~dd903/assets/papers/vldb16-setjoin.pdf)) and overlap set join ([Deng et al. SIGMOD 2018](https://people.cs.rutgers.edu/~dd903/assets/papers/sigmod18.pdf)) algorithms to generate the candidate pool within the time limit. Our approach doesn't rely on heavy data preprocessing or fancy manual rules.

- Dong Deng, Yufei Tao, Guoliang Li: Overlap Set Similarity Joins with Theoretical Guarantees. SIGMOD Conference 2018: 905-920 [bibtex](https://dblp.org/rec/conf/sigmod/DengT018.html?view=bibtex)
- Dong Deng, Guoliang Li, He Wen, Jianhua Feng: An Efficient Partition Based Method for Exact Set Similarity Joins. Proc. VLDB Endow. 9(4): 360-371 (2015) [bibtex](https://dblp.org/rec/journals/pvldb/DengLWF15.html?view=bibtex)

## Execution

```bash
bash run.sh
```
This command will generate output.csv, trace the workflow and generate submission.rpz.

```bash
bash run_without_reprozip.sh
```
This command will generate output.csv directly, without using repro.

## Other Scripts

The `run.sh` and `run_without_reprozip.sh` execute our preprocessing and blocking scripts. 

Specifically, we have two atomic scripts:

    `sed_clean.sh` is our preprocessing script to clean the raw data using sed. It takes X1.csv and X2.csv as input and outputs clean_X1.csv and clean_X2.csv.

    `run_blocking.sh` is our blocking script to execute blocking on the two datasets clean_X1.csv and clean_X2.csv, and generate the final output.csv.

## Executable Program

Our core executable program `blocking` takes 8 parameters to work.

```bash
./blocking SOURCE_CSV OUTPUT_CSV RESULT_NUM GROUP_THRESHOLD JACCARD_THRESHOLD OVERLAP_THRESHOLD POOL_SIZE_FACTOR
```

    `SOURCE_CSV`: string, source dataset for blocking. It must be a CSV file.
    `OUTPUT_CSV`: string, output file name containing blocking results. Each line is a pair of record IDs.
    `RESULT_NUM`:  int, the number of possible duplicate record pairs to generate after blocking
    `GROUP_THRESHOLD`: float in (0,1], jaccard threshold to group different entries
    `JACCARD_THRESHOLD`: float in  (0,1], jaccard threshold to filter out non-candidate pairs
    `OVERLAP_THRESHOLD`: int in [1, inf), number of overlap size to filter out non-candidate pairs
    `POOL_SIZE_FACTOR`: int, determines the size of the candidate pool (= POOL_SIZE_FACTOR * RESULT_NUM), used to avoid memory overflow

An example complete command is: 

`./blocking clean_X1.csv output1.csv 1000000 0.8 0.57 9 15`

