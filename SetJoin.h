
#ifndef _SETJOIN_H_
#define _SETJOIN_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <algorithm>
#include <string.h>
#include <queue>
#include <inttypes.h>
#include <sys/time.h>
#include <cmath>
#include <cstdio>

using namespace std;

#define google_unordered_map unordered_map
#define google_unordered_set unordered_set

#define PACK(x, y) ((x << 32) + y)

#define PRIME 2017
#define EPS 1e-5
#define NEG -1
#define INF 100000000
#define MAX_LINE_LENGTH 1000000
#define CACHE_SIZE 5

extern vector<pair<int, int>> cacheVec;
extern vector<vector<pair<int, int>>> indexVecs;


class SimPairx
{
public:
    int id1;
    int id2;
    double sim;
    SimPairx() {}
    SimPairx(int i1, int i2, double s)
    {
        id1 = i1;
        id2 = i2;
        sim = s;
    }
    
    bool operator<(const SimPairx& rhs) const
    {
        return this->sim > rhs.sim;
    }
};

class SetJoin 
{
public:
    int same, diff;
    double det;
    uint64_t resultNum = 0;
    uint64_t candidateNum = 0;
    uint64_t lengthSum = 0;
    uint64_t listlens = 0;
    int prime_exp[MAX_LINE_LENGTH];
    vector<vector<int>> dataset;
    vector<double> &wordwt;
    vector<double> &recordwt;
    int maxlimit;
    bool has_limit;
    vector<pair<int, int>> result_pairs;

    priority_queue<SimPairx> result_pairs_;

    void get_results()
    {

                while(!result_pairs_.empty())
                {
                    result_pairs.emplace_back(result_pairs_.top().id1, result_pairs_.top().id2);
                    result_pairs_.pop();
                }
    }

    struct invertedList {
        int vec_no, cnt;
        pair<int, int> cache[CACHE_SIZE];

        vector<pair<int, int>>& getVector() const {
            if (cnt <= CACHE_SIZE) {
                cacheVec.assign(cache, cache + cnt);
                return cacheVec;
            } else
                return indexVecs[vec_no];
        }

        void add(pair<int, int> data) {
            if (cnt < CACHE_SIZE) cache[cnt++] = data;
            else {
                if (CACHE_SIZE == cnt) {
                    indexVecs.push_back(vector<pair<int, int>>());
                    vec_no = indexVecs.size() - 1;
                    indexVecs[vec_no].assign(cache, cache + CACHE_SIZE);
                }
                ++cnt;
                indexVecs[vec_no].push_back(data);
            }
        }
    };

    struct invIndexStruct {
        int list_no;
        int* oneList;
        //vector<int> oneList;

        invIndexStruct() : list_no(0) {}
    };

    vector<invertedList> indexLists;
    
    SetJoin(vector<vector<int>> &sorted_records, vector<double> &ww, vector<double> &rw, int ml, bool islimit)
        : wordwt(ww), recordwt(rw)
    {
        dataset = sorted_records;
        maxlimit = ml;
        has_limit = islimit;

        indexVecs.clear();
        cacheVec.clear();
        cacheVec.resize(CACHE_SIZE);

    }

    ~SetJoin()
    {
        cacheVec.clear();
        indexVecs.clear();
    }

    bool overlap(int x, int y, int posx = 0, int posy = 0, int current_overlap = 0);
    void setjoin(double threshold); 
};

#endif