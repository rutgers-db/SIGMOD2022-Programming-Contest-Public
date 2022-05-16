#ifndef _OVLPJOIN_H_
#define _OVLPJOIN_H_

#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <string.h>
#include <iostream>
#include <algorithm>
#include <inttypes.h>
#include <queue>
#include <sys/time.h>
#include <time.h>
#include <unordered_map>
#include <unordered_set>

using namespace std;

#define RATIO 0.005
#define TIMES 200



class SimPairy
{
public:
    int id1;
    int id2;
    double sim;
    SimPairy() {}
    SimPairy(int i1, int i2, double s)
    {
        id1 = i1;
        id2 = i2;
        sim = s;
    }
    
    bool operator<(const SimPairy& rhs) const
    {
        return this->sim > rhs.sim;
    }
};


struct combination;

extern int c;
extern vector<vector<int>> dataset;
extern vector<combination> combs;
extern bool comp_int(const int a, const int b);
extern bool comp_comb(const int a, const int b);
extern bool comp_pair(const pair<int, int> &p1, const int val);

struct combination
{
  int N;
  int id;
  bool completed;
  vector<int> curr;

  combination(int d, int beg)
      : completed(false), id(d), N(dataset[d].size())
  {
    if (N < 1 || c > N)
      completed = true;
    for (auto i = 0; i < c; ++i)
      curr.push_back(beg + 1 + i);
  }

  // compute next combination
  void next()
  {
    int i = c - 1;
    while (i >= 0 && curr[i] == N - c + i)
      --i;
    if (i < 0)
      completed = true;
    else
    {
      int temp = curr[i];
      for (int j = i; j <= c - 1; j++)
        curr[j] = temp + 1 + j - i;
    }
  }

  void print() const
  {
    cout << "combination from " << id << " : ";
    for (auto j = 0; j < c; j++)
      cout << dataset[id][curr[j]] << " ";
    cout << " ----> ";
    for (auto j = 0; j < c; j++)
      cout << curr[j] << " ";
    cout << endl;
  }

  bool stepback(int i)
  {
    if (i == 0)
      return true;
    curr[i - 1]++;
    if (curr[i - 1] + c - 1 - i + 1 >= N)
      return stepback(i - 1);
    for (int j = i; j < c; j++)
      curr[j] = curr[i - 1] + j - i + 1;
    return false;
  }

  void binary(const combination &value)
  {
    auto it = dataset[id].begin() + curr[0];
    for (int i = 0; i < c; i++)
    {
      // find the first one not larger than the value
      it = lower_bound(it, dataset[id].end(), dataset[value.id][value.curr[i]], comp_int);
      // if get the end, we will increase the last one by 1 and set the rest as max
      if (it == dataset[id].end())
      {
        completed = stepback(i);
        return;
        // if we get the same value, we fill in it
      }
      else if (*it == dataset[value.id][value.curr[i]])
      {
        curr[i] = distance(dataset[id].begin(), it);
        // if we get the smaller value, we set the rest as max
      }
      else
      {
        curr[i] = distance(dataset[id].begin(), it);
        if (curr[i] + c - 1 - i >= N)
        {
          completed = stepback(i);
          return;
        }
        for (int j = i + 1; j < c; j++)
          curr[j] = curr[i] + j - i;
        return;
      }
    }
    return;
  }
};

class OvlpJoin 
{
  

public:

  int n;
  int total_eles;
  int alive_id = 0;
  uint64_t heap_op = 0;
  int64_t large_cost = 0;
  int64_t large_est_cost = 0;
  vector<int> heap;
  vector<pair<int, int>> buck;
  vector<vector<int>> records;
  vector<pair<int, int>> idmap;
  unordered_set<int> random_ids;
  vector<vector<pair<int, int>>> ele_lists;
  vector<pair<int, int>> result_pairs;
    priority_queue<SimPairy> result_pairs_;
    vector<double> &wordwt;
    vector<double> &recordwt;
    int maxlimit;
    bool has_limit;



    void get_results()
    {

                while(!result_pairs_.empty())
                {
                    result_pairs.emplace_back(result_pairs_.top().id1, result_pairs_.top().id2);
                    result_pairs_.pop();
                }
    }
    
  double heap_cost;
  double binary_cost;
  int list_max;
  int list_min;
  int64_t list_cost;
  int64_t list_sum;
  int64_t list_sample_num;
  int64_t result_num;
  int64_t candidate_num;

  void overlapjoin(int overlap_threshold);
  bool is_equal(const combination & c1, const combination & c2);
  int estimate();
  int64_t small_estimate(int L, int R);
  int64_t large_estimate(int L, int R);
  void small_case(int L, int R);
  void large_case(int L, int R);

  uint64_t getListCost();
  int divide(int nL);

  OvlpJoin(vector<vector<int>> &sorted_records, vector<double> &ww, vector<double> &rw, int ml, bool islimit)
        : wordwt(ww), recordwt(rw)
    {
        maxlimit = ml;
        has_limit = islimit;

    // reset everything
    c = 0;
    n = 0;
    total_eles =0;
    alive_id = 0;
    heap_op = 0;
    large_cost = 0;
    large_est_cost = 0;
    
    heap.clear();
    buck.clear();
    dataset.clear();
    records.clear();
    idmap.clear();
    random_ids.clear();
    ele_lists.clear();
    result_pairs.clear();

    combs.clear();

    heap_cost = 0;
    binary_cost = 0;
    list_max = 0;
    list_min = 0;
    list_cost = 0;
    list_sum = 0;
    list_sample_num = 0;
    result_num = 0;
    candidate_num = 0;
    
    records = sorted_records;
  }
};
#endif
