
#include "OvlpJoin.h"


double similarityx(vector<int> &v1, vector<int> &v2, double w1, double w2, const vector<double> &ww)
{
  double ovlp = 0.0; 
  auto it1 = v1.begin();
  auto it2 = v2.begin();

  while (it1 != v1.end() && it2 != v2.end()) {
    if (*it1 == *it2) {
      ovlp += ww[*it1];
      ++it1, ++it2;
    } else {
      if (*it1 < *it2) ++it1;
      else ++it2;
    }
  }
  // jaccard
  double jac = ovlp / (w1 + w2 - ovlp);
  double olp = ovlp / min(w1, w2);
  
  // return jac + olp;
  // rank by L2-norm?
  // return jac * jac + olp;
  return olp;
}


double similarityy(vector<int> &v1, vector<int> &v2, double w1, double w2, const vector<double> &ww, double t1, double t2) 
{
  int common = 0;
  double ovlp = 0.0; 
  auto it1 = v1.begin();
  auto it2 = v2.begin();

  while (it1 != v1.end() && it2 != v2.end()) {
    if (*it1 == *it2) {
      ovlp += ww[*it1];
      ++common;
      ++it1, ++it2;
    } else {
      if (*it1 < *it2) ++it1;
      else ++it2;
    }   
  }
  
  // rank by MAX
  double sim = 0;
  // jaccard
  if (ovlp / (w1 + w2 - ovlp) >= t1) 
    sim = sim + ovlp / (w1 + w2 - ovlp);

  // overlap
  if (common >= t2) 
    sim = sim + ovlp / min(w1, w2);

  return sim; 
}

int c;
vector<vector<int>> dataset;
vector<combination> combs;

int64_t nchoosek(int64_t n, int64_t k) {
  if (k == 0) return 1;
  return (n * nchoosek(n - 1, k - 1)) / k;
}

void OvlpJoin::small_case(int L, int R) {
  if (L >= R) return;
  --c;

  timeval beg, mid, mid1, end;
  gettimeofday(&beg, NULL);

  cout << " number of small sets: " << R - L << endl;

  list_cost = 0;
  heap_cost = 0;
  binary_cost = 0;
  vector<vector<int>> res_lists;

  gettimeofday(&mid, NULL);

  for (auto idx = total_eles - 1; idx >= 0; idx--) {
    if (ele_lists[idx].size() < 2) continue;
    vector<pair<int,int>> & vec = ele_lists[idx];
    int size = distance(vec.begin(), lower_bound(vec.begin(), vec.end(), L, comp_pair));
    if (vec.size() <= size + 1) continue;

    // initialize heap and combinations
    heap.clear();
    combs.clear();
    int heap_size = 0;
    for (auto i = size; i < vec.size(); i++) {
      if ((int)(dataset[vec[i].first].size()) - 1 - vec[i].second < c) continue;
      heap.push_back(heap_size++);
      combs.push_back(combination(vec[i].first, vec[i].second));
    }

    if (heap_size < 2) continue;
    make_heap(heap.begin(), heap.end(), comp_comb);
    heap_cost += (3 * c * heap_size);
    // cout << heap_size << " initial: " << heap_cost << endl;

    // pop heaps
    vector<int> inv_list;
    while (heap_size > 1) {
      inv_list.clear();
      do {
        heap_cost += (c * log2(heap_size) + c);
        // cout << heap_size << " " << heap_cost << endl;
        pop_heap(heap.begin(), heap.begin() + heap_size, comp_comb);
        --heap_size;
        inv_list.push_back(combs[heap[heap_size]].id);
      } while (heap_size > 0 && is_equal(combs[heap[heap_size]], combs[heap.front()]));

      if (inv_list.size() > 1) {
        list_cost += ((inv_list.size() - 1) * (int64_t)inv_list.size() / 2);
        res_lists.push_back(std::move(inv_list));
      }

      if (heap_size == 0) break;

      for (auto i = heap_size; i < heap.size(); ++i) {
        combs[heap[i]].binary(combs[heap.front()]);
        binary_cost += (c * log2(dataset[combs[heap[i]].id].size()));
      }

      int comp_num = 0;
      for (auto i = heap_size; i < heap.size(); ++i) {
        if (combs[heap[i]].completed)
          ++comp_num;
        else if (comp_num > 0)
          heap[i - comp_num] = heap[i];
      }

      for (auto i = heap_size; i < (int)heap.size() - comp_num; i++) {
        push_heap(heap.begin(), heap.begin() + i + 1, comp_comb);
        heap_cost += (c * log2(i));
      }
      while (comp_num-- > 0)
        heap.pop_back();
      heap_size = heap.size();
    }
  }

  cout << "Res lists num: " << res_lists.size() << endl;

  gettimeofday(&mid1, NULL);

  vector<vector<int>> id_lists(n);
  for (auto i = 0; i < res_lists.size(); i++) {
    for (auto j = 0; j < res_lists[i].size(); j++)
      id_lists[res_lists[i][j]].push_back(i);
  }

  vector<int> results(n, -1);
  for (auto i = n - 1; i >= 0; i--) {
    if (id_lists[i].empty()) continue;
    for (auto j = 0; j < id_lists[i].size(); j++) {
      res_lists[id_lists[i][j]].pop_back();
      for (auto k = 0; k < res_lists[id_lists[i][j]].size(); k++) {
        if (results[res_lists[id_lists[i][j]][k]] != i) {
          // cout << idmap[i].first << " " << idmap[res_lists[id_lists[i][j]][k]].first << endl;
          results[res_lists[id_lists[i][j]][k]] = i;
          int idd1 = idmap[i].first; 
          int idd2 = idmap[res_lists[id_lists[i][j]][k]].first;
            if (has_limit)
                        {
                            double sim = similarityx(dataset[idd1], dataset[idd2], recordwt[idd1], recordwt[idd2], wordwt);
                            if (result_pairs_.size() > maxlimit)
                                result_pairs_.pop();
                            result_pairs_.emplace(idd1, idd2, sim);
                        } else {
                            result_pairs.emplace_back(idd1, idd2);
                        }

          ++result_num;
        }
      }
    }
  }
  ++c;
  gettimeofday(&end, NULL);
  cout << " small p1 : " << mid.tv_sec - beg.tv_sec + (mid.tv_usec - beg.tv_usec) / 1e6 << endl;
  cout << " small p2 : " << mid1.tv_sec - mid.tv_sec + (mid1.tv_usec - mid.tv_usec) / 1e6 << endl;
  cout << " small p3 : " << end.tv_sec - mid1.tv_sec + (end.tv_usec - mid1.tv_usec) / 1e6 << endl;
  cout << " heap, binary, list costs : " << heap_cost << " " << binary_cost << " " << list_cost << endl;
}

void OvlpJoin::large_case(int L, int R) {
  timeval beg, mid, end;
  gettimeofday(&beg, NULL);

  vector<vector<int>> ele(total_eles);

  for (int i = n - 1; i >= R; --i)
    for (auto x : dataset[i])
      ele[x].push_back(i);

  gettimeofday(&mid, NULL);
  vector<int> bucket;

  for (int i = R - 1; i >= L; --i) {

    int count = 0;
    for (auto x : dataset[i]) count += ele[x].size();
    large_cost += count;

    if (count > 0.2 * n) {
      // n + count * value + count * if
      bucket.assign(n, 0);
      for (auto x : dataset[i]) {
        for (auto id : ele[x]) {
          if (++bucket[id] == c) {
            // cout << idmap[i].first << " " << idmap[id].first << endl;
              int idd1 = idmap[i].first; 
              int idd2= idmap[id].first;


                        if (has_limit)
                        {
                            double sim = similarityx(dataset[idd1], dataset[idd2], recordwt[idd1], recordwt[idd2], wordwt);
                            if (result_pairs_.size() > maxlimit)
                                result_pairs_.pop();
                            result_pairs_.emplace(idd1, idd2, sim);
                        } else {
                            result_pairs.emplace_back(idd1, idd2);
                        }


            ++result_num;
          }
        }
      }
    } else {
      // count * if + count * value + count * if or value
      alive_id++;
      for (auto x : dataset[i]) {
        for (auto id : ele[x]) {
          if (buck[id].first != alive_id) {
            buck[id].first = alive_id;
            buck[id].second = 1;
          } else {
            if (++buck[id].second == c) {
              // cout << idmap[i].first << " " << idmap[id].first << endl;
              int idd1 = idmap[i].first;
              int idd2 = idmap[id].first;


                        if (has_limit)
                        {
                            double sim = similarityx(dataset[idd1], dataset[idd2], recordwt[idd1], recordwt[idd2], wordwt);
                            if (result_pairs_.size() > maxlimit)
                                result_pairs_.pop();
                            result_pairs_.emplace(idd1, idd2, sim);
                        } else {
                            result_pairs.emplace_back(idd1, idd2);
                        }
              ++result_num;
            }
          }
        }
      }
    }
    for (auto x : dataset[i])
      ele[x].push_back(i);
  }
  gettimeofday(&end, NULL);
  cout << " large p1 : " << mid.tv_sec - beg.tv_sec + (mid.tv_usec - beg.tv_usec) / 1e6 << endl;
  cout << " large p2 : " << end.tv_sec - mid.tv_sec + (end.tv_usec - mid.tv_usec) / 1e6 << endl;
}

// return the bound position
int OvlpJoin::estimate() {

  // get random elements for sampling
  while (random_ids.size() < total_eles * RATIO)
    random_ids.insert(rand() % total_eles);

  int64_t small, large;
  int min_size = dataset.back().size();
  int max_size = dataset.front().size();
  auto bound = (min_size <= c ? c : min_size);
  int pos = divide(bound);
  int prev_pos = pos;
  int64_t prev_large = large_estimate(0, pos);
  int64_t prev_small = small_estimate(pos, n);
  ++bound;

  for (; bound <= max_size; bound++) {

    pos = divide(bound);
    if (pos == prev_pos) continue;
    cout << endl << "size boud: " << bound << endl;
    cout << "larg numb: " << pos << endl;
    cout << "smal numb: " << n - pos << endl;

    large = large_estimate(0, pos);
    small = small_estimate(pos, n);

    cout << "heap cost: " << heap_cost * TIMES << endl;
    cout << "biny cost: " << binary_cost * TIMES << endl;
    cout << "list cost: " << list_cost << endl; 
    cout << "smal cost: " << small << endl;
    cout << "larg cost: " << large << endl;

    if (small - prev_small > 1.2 * (prev_large - large)) return prev_pos;

    prev_pos = pos;
    prev_large = large;
    prev_small = small;
  }
  return prev_pos;
}

void OvlpJoin::overlapjoin(int overlap_threshold)
{
  srand(time(NULL));
  
  timeval starting, ending, s1, t1, s2, t2;
  timeval time1, time2, time3, time4;

  gettimeofday(&starting, NULL);

  c = overlap_threshold; // get threshold
  n = records.size(); // get number of records
  buck.assign(n, make_pair(0, 0)); // for counting

  vector<pair<int, int>> eles;
  unordered_map<int, vector<int>> ele;

  for (int i = 0; i < records.size(); i++) {
    if (records[i].size() < c) continue;  // remove records with size smaller than c
    for (int j = 0; j < records[i].size(); j++)  // build inverted index
      ele[records[i][j]].push_back(i);
  }

  for (auto it = ele.begin(); it != ele.end(); it++)
    eles.push_back(make_pair(it->first, it->second.size()));  // build element frequency table

  // get global order: frequency increasing order
  sort(eles.begin(), eles.end(), [](const pair<int, int> &p1, const pair<int, int> &p2) {
    return p1.second < p2.second;
  });

  // container initialize
  dataset.resize(n);

  // sort elements by its global order: frequence increasing order
  // remove widow word
  // encode elements in decreasing order
  total_eles = eles.size();
  for (auto i = 0; i < int(eles.size()); ++i) {
    if (eles[i].second < 2) continue;
    for (auto j = ele[eles[i].first].begin(); j != ele[eles[i].first].end(); j++)
      dataset[*j].push_back(total_eles - i - 1);
  }

  gettimeofday(&time1, NULL);
  cout << "Initial Time: " << time1.tv_sec - starting.tv_sec + (time1.tv_usec - starting.tv_usec) / 1e6 << endl;

  // ****** cost model for prefix length selection ******
  // remove short records
  for (auto i = 0; i < n; i++)
    if (dataset[i].size() < c) dataset[i].clear();

  // create id mappings: from sorted to origin
  for (auto i = 0; i < n; i++)
    idmap.push_back(make_pair(i, dataset[i].size()));

  // sort records by length in decreasing order
  sort(idmap.begin(), idmap.end(), [] (const pair<int, int> & a, const pair<int, int> & b) {
      return a.second > b.second;
  });
  sort(dataset.begin(), dataset.end(), [] (const vector<int>& a, const vector<int>& b) {
      return a.size() > b.size();
  });
  cout << " largest set: " << dataset.front().size() << " smallest set: " << dataset.back().size() << endl;

  // build real inverted index
  ele_lists.resize(total_eles);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < dataset[i].size(); j++)
      ele_lists[dataset[i][j]].push_back(make_pair(i,j));

  gettimeofday(&time3, NULL);
  cout << "Transform Time: " << time3.tv_sec - time1.tv_sec + (time3.tv_usec - time1.tv_usec) / 1e6 << endl;

  // ****** cost model for boundary selection ******
  int nL = estimate();
  int nP = nL;
  cout << " large sets: " << nP << " small sets: " << n - nP << endl;

  gettimeofday(&time4, NULL);
  cout << "Estimation Time: " << time4.tv_sec - time3.tv_sec + (time4.tv_usec - time3.tv_usec) / 1e6 << endl;

  // ****** conduct joining ******  
  result_num = 0;
  candidate_num= 0;

  gettimeofday(&s1, NULL);
  large_case(0, nP);
  gettimeofday(&t1, NULL);

  gettimeofday(&s2, NULL);
  small_case(nP, n);
  gettimeofday(&t2, NULL);

  gettimeofday(&ending, NULL);
  cout << "Join Time: " << ending.tv_sec - time4.tv_sec + (ending.tv_usec - time4.tv_usec) / 1e6 << endl;
  cout << "  large Time: " << t1.tv_sec - s1.tv_sec + (t1.tv_usec - s1.tv_usec) / 1e6 << endl;
  cout << "  small Time: " << t2.tv_sec - s2.tv_sec + (t2.tv_usec - s2.tv_usec) / 1e6 << endl;
  cout << "All Time: " << ending.tv_sec - starting.tv_sec + (ending.tv_usec - starting.tv_usec) / 1e6 << endl;
  cout << "Result Num: " << result_num << endl;
  cout << "  large cost: " << large_cost << " small cost: " << heap_cost + list_cost + binary_cost << endl;
}

int64_t OvlpJoin::small_estimate(int L, int R) {
  if (L >= R) return 0;
  
  timeval beg, mid, mid1, end;
  gettimeofday(&beg, NULL);

  int total_num = R - L;
  int sample_time = (R - L);
  double ratio =  (total_num - 1) * 1.0 / sample_time * total_num / 2;
  // cout << "sample ratio: " << ratio << endl;
  int r1, r2;
  int64_t pair_num = 0;
  for (auto i = 0; i < sample_time; i++) {
    do {
      r1 = rand() % (R - L) + L;
      r2 = rand() % (R - L) + L;
    } while (r1 == r2);
    int start1 = 0;
    int start2 = 0;
    int overlap = 0;
    while (start1 < dataset[r1].size() && start2 < dataset[r2].size()) {
      if (dataset[r1][start1] == dataset[r2][start2]) {
        ++start1, ++start2;
        overlap++;
      } else {
        if (dataset[r1][start1] > dataset[r2][start2]) ++start1;
        else ++start2;
      }
    }
    if (overlap >= c){
      // cout << overlap << " " << c << endl;
      pair_num += nchoosek(overlap, c);
      // cout << pair_num << endl;
    }
  }
  list_cost = pair_num * ratio;

  --c;

  heap_cost = 0;
  binary_cost = 0;

  gettimeofday(&mid, NULL);

  for (auto sit = random_ids.begin(); sit != random_ids.end(); ++sit) {
    auto idx = *sit;

    vector<pair<int,int>> & vec = ele_lists[idx];
    int size = distance(vec.begin(), lower_bound(vec.begin(), vec.end(), L, comp_pair));
    if (vec.size() <= size + 1) continue;

    heap.clear();
    combs.clear();
    int heap_size = 0;
    for (auto i = size; i < vec.size(); i++) {
      if ((int)(dataset[vec[i].first].size()) - 1 - vec[i].second < c) continue;
      heap.push_back(heap_size++);
      combs.push_back(combination(vec[i].first, vec[i].second));
    }

    if (heap_size < 2) continue;

    make_heap(heap.begin(), heap.end(), comp_comb);
    heap_cost += (3 * c * heap_size);

    while (heap_size > 1) {
      do {
        ++heap_op;
        heap_cost += (c * log2(heap_size) + c);
        pop_heap(heap.begin(), heap.begin() + heap_size, comp_comb);
        --heap_size;
      } while (heap_size > 0 && is_equal(combs[heap[heap_size]], combs[heap.front()]));

      if (heap_size == 0) break;

      for (auto i = heap_size; i < heap.size(); ++i) {
        combs[heap[i]].binary(combs[heap.front()]);
        binary_cost += (c * log2(dataset[combs[heap[i]].id].size()));
      }

      int comp_num = 0;
      for (auto i = heap_size; i < heap.size(); ++i) {
        if (combs[heap[i]].completed)
          ++comp_num;
        else if (comp_num > 0)
          heap[i - comp_num] = heap[i];
      }

      for (auto i = heap_size; i < (int)heap.size() - comp_num; i++) {
        push_heap(heap.begin(), heap.begin() + i + 1, comp_comb);
        heap_cost += (c * log2(i));
      }
      while (comp_num-- > 0)
        heap.pop_back();
      heap_size = heap.size();
    }
  }

  ++c;

  gettimeofday(&end, NULL);
  cout << " small est time p1 : " << mid.tv_sec - beg.tv_sec + (mid.tv_usec - beg.tv_usec) / 1e6 << endl;
  cout << " small est time p2 : " << end.tv_sec - mid.tv_sec + (end.tv_usec - mid.tv_usec) / 1e6 << endl;
  return binary_cost * TIMES + heap_cost * TIMES + list_cost;
}

/*
int64_t large_estimate(int L, int R) {
  int64_t ret = 0;
  for (int x = 0; x < ele_lists.size(); x++) {
    int large_num = distance(ele_lists[x].begin(), lower_bound(ele_lists[x].begin(), ele_lists[x].end(), R, comp_pair));
    int all_num = ele_lists[x].size();
    ret += (large_num * all_num - large_num * (large_num - 1) / 2);
  }
  large_est_cost = ret;
  return ret;
}
*/

int64_t OvlpJoin::large_estimate(int L, int R) {
  timeval beg, end;
  gettimeofday(&beg, NULL);
  vector<int> count(total_eles);
  for (int i = n - 1; i >= R; --i)
    for (auto x : dataset[i])
      ++count[x];

  int64_t ret = 0;
  for (int i = R - 1; i >= L; --i) {
    for (auto x : dataset[i]) {
        ++count[x];
        ret += count[x];
    }
  }
  large_est_cost = ret;
  gettimeofday(&end, NULL);
  cout << " large est time : " << end.tv_sec - beg.tv_sec + (end.tv_usec - beg.tv_usec) / 1e6 << endl;
  return ret;
}

uint64_t OvlpJoin::getListCost() {
  return (list_cost - list_sum) / 2 * TIMES;
  // return (list_cost - (list_sum * 1.0 / list_sample_num) * (list_sum * 1.0 / list_sample_num) - list_sum) / 2 * TIMES;
  // return (list_cost * list_sample_num * 1.0 / (list_sample_num - 1) - list_sum * 1.0 / list_sample_num * list_sum / (list_sample_num - 1) - list_sum) / 2 * TIMES; 
}

// find first set with size smaller or equal to nL
int OvlpJoin::divide(int nL) {
  int l = 0, r = n;
  while (l < r) {
    int m = (l + r) >> 1;
    if (dataset[m].size() > nL) l = m + 1;
    else r = m;
  }
  return r;
}

bool comp_pair(const pair<int, int> &p1, const int val) {
  return p1.first < val;
}

bool comp_int(const int a, const int b) {
  return a > b;
}

bool comp_comb(const int a, const int b) {
  combination & c1 = combs[a];
  combination & c2 = combs[b];
  for (int i = 0; i < c; i++) {
    if (dataset[c1.id][c1.curr[i]] > dataset[c2.id][c2.curr[i]])
      return false;
    else if (dataset[c1.id][c1.curr[i]] < dataset[c2.id][c2.curr[i]])
      return true;
  }
  return c1.id > c2.id;
}

bool OvlpJoin::is_equal(const combination & c1, const combination & c2) {
  for (int i = 0; i < c; i++) {
    if (dataset[c1.id][c1.curr[i]] != dataset[c2.id][c2.curr[i]])
      return false;
  }
  return true;
}
