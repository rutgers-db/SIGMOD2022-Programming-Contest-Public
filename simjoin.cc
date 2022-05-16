#include <iostream>
#include <fstream>
#include <vector> 
#include <string>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <math.h>
#include <numeric>
#include "CSVReader.h"
#include "Table.h"
#include "SetJoin.h"
#include "OvlpJoin.h"

int HEAPRATE;

#ifdef IN_PARALLEL
#include <omp.h>
#endif

class DSU {
public:
    DSU(int n) {
        data.resize(n);
        iota(begin(data), end(data), 0);        
        sizes.resize(n);
        std::fill(begin(sizes), end(sizes), 1);
        // pair_num = 0;
    }
    
    // path halving
    int find(int x) {
        while (x != data[x]) {
            data[x] = data[data[x]];
            x = data[x];
        }
        return x;
    }
    
    bool join(int x, int y) {
        int px = find(x);
        int py = find(y);
        if (px == py) return false;
        // pair_num += (sizes[py] * sizes[px]);

        if (sizes[px] > sizes[py]) 
        {
            data[py] = px;    
            sizes[px] += sizes[py];
        }
        else 
        {
            data[px] = py;    
            sizes[py] += sizes[px];
        }
        return true;
    }
    

    /*
    void reset_pair_num() {
        pair_num = 0;
    }

    int get_pair_num() {
        return pair_num;
    }

    bool over_limit(int limit) {
        return pair_num >= limit;
    }
    */

private:
    // int pair_num = 0;
    vector<int> data;
    vector<int> sizes;
};

void strToTokens(const string &s, vector<string> &res, const string &delims) {
	string::size_type begIdx, endIdx;
	begIdx = s.find_first_not_of(delims);
	while (begIdx != string::npos) {
		endIdx = s.find_first_of(delims, begIdx);
		if (endIdx == string::npos)
			endIdx = s.length();
		res.push_back(s.substr(begIdx, endIdx - begIdx));
		begIdx = s.find_first_not_of(delims, endIdx);
	}
}

class SimPair
{
public:
    int id1;
    int id2;
    double sim;
    SimPair() {
        id1 = 0;
        id2 = 0;
        sim = 0.0;
    }
    SimPair(int i1, int i2, double s)
    {
        id1 = i1;
        id2 = i2;
        sim = s;
    }
};

double jaccard(vector<int> &v1, vector<int> &v2) 
{
  int ovlp = 0.0; 
  auto it1 = v1.begin();
  auto it2 = v2.begin();

  while (it1 != v1.end() && it2 != v2.end()) {
    if (*it1 == *it2) {
      ovlp++;
      ++it1, ++it2;
    } else {
      if (*it1 < *it2) ++it1;
      else ++it2;
    }
  }

  return ovlp * 1.0 / (v1.size() + v2.size() - ovlp);
}


double overlap(vector<int> &v1, vector<int> &v2) 
{
  int ovlp = 0.0; 
  auto it1 = v1.begin();
  auto it2 = v2.begin();

  while (it1 != v1.end() && it2 != v2.end()) {
    if (*it1 == *it2) {
      ovlp++;
      ++it1, ++it2;
    } else {
      if (*it1 < *it2) ++it1;
      else ++it2;
    }
  }

  return ovlp * 1.0 / min(v1.size(), v2.size());
}

double similarity1(vector<int> &v1, vector<int> &v2, double w1, double w2, const vector<double> &ww)
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
  return jac + olp;
}


double similarity2(vector<int> &v1, vector<int> &v2, double w1, double w2, const vector<double> &ww, double t1, double t2) 
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

void txtnorm(Table &table, vector<int> &record_ids, int &word_num, vector<double> &weight, 
    vector<int> &idmap, vector<vector<int>> &sorted_records, vector<string> &all_brands, vector<int> &record_brands)
{
    int rec_num = table.rows.size();

    // get description    
    vector<vector<string>> bows(rec_num);
    string delim = " ,-_;.?!|()\\\t\n\r";
    record_ids.clear();
    record_ids.resize(table.rows.size());

#ifdef IN_PARALLEL
    #pragma omp parallel for
#endif
    for (auto rid = 0; rid < table.rows.size(); rid++)
    {
        auto &row = table.rows[rid];

        string &row_id = row[0];
        record_ids[rid] = atoi(row_id.c_str());

        auto &str = row[1];
        for (int i = 0; i < str.length(); i++)
            str[i] = tolower(str[i]);
        
        auto &bag_of_words = bows[rid];
        strToTokens(str, bag_of_words, delim);

        // sort
        sort(bag_of_words.begin(), bag_of_words.end());

        // unique
        auto uniq_it = unique(bag_of_words.begin(), bag_of_words.end()); 
        bag_of_words.resize(distance(bag_of_words.begin(), uniq_it)); // 10 20 30 20 10
    }

    // get all the words; single thread
    unordered_map<string, vector<int>> inv_index;
    for (auto rid = 0; rid < bows.size(); rid++)
        for (auto &word_str : bows[rid])
            inv_index[word_str].emplace_back(rid);
    word_num = inv_index.size();

    // get brand info
    for (int brand_id = 0; brand_id < all_brands.size(); brand_id++)
    {
        if (inv_index.find(all_brands[brand_id]) != inv_index.end())
            for (auto rid : inv_index[all_brands[brand_id]])
                record_brands[rid] = brand_id;
    }

    // single thread
    vector<pair<int, string>> tokens;
    for (auto &entry : inv_index)
        tokens.emplace_back(entry.second.size(), entry.first);

    // in ascending order of frequency
    sort(tokens.begin(), tokens.end(), [](const pair<int, string> &p1, const pair<int, string> &p2){
        return p1.first < p2.first;
    });

    // print info
    /*
    cout << "dictionary (should be increasing frequnency order) " << endl;
    for (auto &entry : tokens)
        cout << entry.first << " " << entry.second << endl;
    cout << " ===== end of dict===" << endl;
    */

    vector<vector<int>> dataset(rec_num);
    // single thread
    for (int wid = 0; wid < tokens.size(); wid++)
    {
        auto &word = tokens[wid].second;
        for (auto rid : inv_index[word])
            dataset[rid].emplace_back(wid);
        
        int freq = tokens[wid].first;
        weight.emplace_back(log10(rec_num * 1.0 / freq));
    }

    idmap.clear();
    for (auto i = 0; i < rec_num; i++)
        idmap.emplace_back(i);
    
    // sort dataset by size first, element second, id third
    sort(idmap.begin(), idmap.end(), [&dataset](const int &id1, const int &id2) {
        int dsize1 = dataset[id1].size();
        int dsize2 = dataset[id2].size();
        if (dsize1 < dsize2)
            return true;
        else if (dsize1 > dsize2)
            return false;
        else {
            for (int i = 0; i < dsize1; i++)
            {
                if (dataset[id1][i] < dataset[id2][i])
                    return true;
                else if (dataset[id1][i] > dataset[id2][i])  
                    return false;
            }
            if (id1 < id2)
                return true;
            else
                return false;
        }
    });

    /********** w/o deduplication *********/
    for (int i = 0; i < idmap.size(); i++)
        sorted_records.emplace_back(dataset[idmap[i]]);
    /********** w/o deduplication **********/
    
    // output statistics
    int total_words = 0;
    for (int i = 0; i < rec_num; i++)
        total_words += sorted_records[i].size();
    cout << "the average number of tokens in each record is: " << total_words * 1.0 / rec_num << endl;
    cout << "the shortest length is " << sorted_records.front().size() << endl;
    cout << "the longest length is " << sorted_records.back().size() << endl;
}


void clustering(double clustering_threshold, vector<vector<int>> &sorted_records,
    vector<vector<int>> &cidmap, vector<vector<int>> &returned_records)
{
    vector<SimPair> all_simpairs;
    vector<double> dummy1;
    vector<double> dummy2;
    SetJoin joiner(sorted_records, dummy1, dummy2, 0, false);

    joiner.setjoin(clustering_threshold);

    int rec_num = sorted_records.size();

    DSU dsu(rec_num);
    for (auto &entry : joiner.result_pairs)
    {        
        int id1 = entry.first;
        int id2 = entry.second;
        dsu.join(id1, id2);
    }

    unordered_map<int, vector<int>> clusters;
    vector<int> tmp_distinct_keys;
    for (int id = 0; id < rec_num; id++)
    {
        int cid = dsu.find(id);

        if (clusters.find(cid) == clusters.end())
            tmp_distinct_keys.emplace_back(cid);

        clusters[cid].emplace_back(id);
    }

    // merge records
    vector<vector<int>> sids(clusters.size());
    vector<vector<int>> merged_records(clusters.size());

#ifdef IN_PARALLEL
    #pragma omp parallel for
#endif
    for (int cnt = 0; cnt < tmp_distinct_keys.size(); cnt++)
    {
        auto &one_cluster = clusters[tmp_distinct_keys[cnt]];
        auto &new_record = merged_records[cnt];
        

        for (auto sid : one_cluster)
        {
            auto &record = sorted_records[sid];
            new_record.insert(new_record.end(), record.begin(), record.end());
            sids[cnt].emplace_back(sid);
        }

        if (one_cluster.size() > 1)
        {
            // sort
            sort(new_record.begin(), new_record.end());
            // unique
            auto uniq_it = unique(new_record.begin(), new_record.end()); 
            new_record.resize(distance(new_record.begin(), uniq_it)); // 10 20 30 20 10   
        }
    }

    vector<int> tmp_idmap;
    for (auto i = 0; i < merged_records.size(); i++)
        tmp_idmap.emplace_back(i);
    
    // sort dataset by size first, element second, id third
    sort(tmp_idmap.begin(), tmp_idmap.end(), [&merged_records](const int &id1, const int &id2) {
        int dsize1 = merged_records[id1].size();
        int dsize2 = merged_records[id2].size();
        if (dsize1 < dsize2)
            return true;
        else if (dsize1 > dsize2)
            return false;
        else {
            for (int i = 0; i < dsize1; i++)
            {
                if (merged_records[id1][i] < merged_records[id2][i])
                    return true;
                else if (merged_records[id1][i] > merged_records[id2][i])  
                    return false;
            }
            if (id1 < id2)
                return true;
            else
                return false;
        }
    });


    /********** w/o deduplication *********/
    returned_records.clear();
    returned_records.reserve(merged_records.size());
    cidmap.clear();
    cidmap.reserve(sids.size());
    for (int i = 0; i < tmp_idmap.size(); i++)
    {
        int tmp_id = tmp_idmap[i];
        returned_records.emplace_back(merged_records[tmp_id]);
        cidmap.emplace_back(sids[tmp_id]);
    }
    /********** w/o deduplication **********/
}

void postprocessing(vector<int> &record_ids, vector<double> &weight, vector<int> &idmap, vector<vector<int>> &cidmap,
    vector<vector<int>> &sorted_records, vector<vector<int>> &merged_records, int limits, ofstream &output, 
    bool use_jaccard, double jaccard_threshold, bool use_overlap, double overlap_threshold, const vector<int> &record_brands)
{
    vector<double> record_wt(sorted_records.size(), 0.0);
#ifdef IN_PARALLEL
    #pragma omp parallel for
#endif 
    for (int i = 0; i < merged_records.size(); i++)
        for (auto wid : merged_records[i])
            record_wt[i] += weight[wid];

#ifdef IN_PARALLEL
    vector<pair<int, int>> jaccard_simpairs;
    vector<pair<int, int>> overlap_simpairs;    
    #pragma omp parallel sections
    {
        #pragma omp section
        {   
            printf ("id = %d, \n", omp_get_thread_num());
            if (use_jaccard)
            {
                //  JACCARD SIMILARITY JOIN
                SetJoin joiner(merged_records, weight, record_wt, HEAPRATE * limits, true);
                joiner.setjoin(jaccard_threshold);
                joiner.get_results();
                for (auto &entry : joiner.result_pairs)
                {
                    int id1 = entry.first;
                    int id2 = entry.second;
                    if (id1 < id2)
                        jaccard_simpairs.emplace_back(id1, id2);
                    else
                        jaccard_simpairs.emplace_back(id2, id1);
                }
            }
        }
        #pragma omp section
        { 
            printf ("id = %d, \n", omp_get_thread_num());
            if (use_overlap)
            {
                // OVERLAP SIMILARITY JOIN
                OvlpJoin joiner(merged_records, weight, record_wt, HEAPRATE * limits, true);
                joiner.overlapjoin(overlap_threshold);
                joiner.get_results();
                
                for (auto &entry : joiner.result_pairs)
                {
                    int id1 = entry.first;
                    int id2 = entry.second;
                    if (id1 < id2)
                        overlap_simpairs.emplace_back(id1, id2);
                    else
                        overlap_simpairs.emplace_back(id2, id1);
                }
            }
        }
    }

    vector<pair<int, int>> &all_simpairs = jaccard_simpairs;
    all_simpairs.insert(all_simpairs.end(), overlap_simpairs.begin(), overlap_simpairs.end());

#else

    vector<pair<int, int>> all_simpairs;
    if (use_jaccard)
    {
        //  JACCARD SIMILARITY JOIN
        SetJoin joiner(merged_records, weight, record_wt, HEAPRATE * limits, true);
        joiner.setjoin(jaccard_threshold);
        joiner.get_results();
        for (auto &entry : joiner.result_pairs)
        {
            int id1 = entry.first;
            int id2 = entry.second;
            if (id1 < id2)
                all_simpairs.emplace_back(id1, id2);
            else
                all_simpairs.emplace_back(id2, id1);
        }
    }

    if (use_overlap)
    {
        // OVERLAP SIMILARITY JOIN
        OvlpJoin joiner(merged_records, weight, record_wt, HEAPRATE * limits, true);
        joiner.overlapjoin(overlap_threshold);
        joiner.get_results();

        for (auto &entry : joiner.result_pairs)
        {
            int id1 = entry.first;
            int id2 = entry.second;
            if (id1 < id2)
                all_simpairs.emplace_back(id1, id2);
            else
                all_simpairs.emplace_back(id2, id1);
        }
    }

#endif 

    sort(all_simpairs.begin(), all_simpairs.end(), [](const pair<int, int> &p1, const pair<int, int> &p2) {
        if (p1.first < p2.first)
            return true;
        else if (p1.first > p2.first)
            return false;
        else if (p1.second < p2.second)
            return true;
        else
            return false;
    });

    // unique
    auto uniq_it = unique(all_simpairs.begin(), all_simpairs.end()); 
    all_simpairs.resize(distance(all_simpairs.begin(), uniq_it)); // 10 20 30 20 10

    /*********** cluster expanding ******************/
    vector<SimPair> simpairs;
    for (auto &entry : cidmap)
        for (auto i = 0; i < entry.size(); i++)
            for (auto j = i + 1; j < entry.size(); j++)
                simpairs.emplace_back(entry[i], entry[j], 0.0);

    for (auto &entry : all_simpairs) 
    {
        for (int sid1 : cidmap[entry.first])
        {
            for (int sid2 : cidmap[entry.second])
            {
                // cout << "#res: " << sid1 << " " << sid2 << endl;
                simpairs.emplace_back(sid1, sid2, 0.0);
            }
        }
    }
    /*********** cluster expanding ******************/
    record_wt.clear();
    record_wt.assign(sorted_records.size(), 0.0);
#ifdef IN_PARALLEL
    #pragma omp parallel for
#endif 
    for (int i = 0; i < sorted_records.size(); i++)
        for (auto wid : sorted_records[i])
            record_wt[i] += weight[wid];

#ifdef IN_PARALLEL
    #pragma omp parallel for
#endif 
    for (auto &entry : simpairs) 
    {
        int id1 = entry.id1;
        int id2 = entry.id2;

        // check brands
        int rowid1 = idmap[id1];
        int rowid2 = idmap[id2];

        if (record_brands[rowid1] >= 0 && record_brands[rowid2] >= 0 && record_brands[rowid1] != record_brands[rowid2])
            entry.sim = EPS;
        else {
            // entry.sim = overlap(sorted_records[id1], sorted_records[id2]);
            // entry.sim += jaccard(sorted_records[id1], sorted_records[id2]);
            //
            if (limits > 1500000)
                entry.sim = similarity2(sorted_records[id1], sorted_records[id2], record_wt[id1], record_wt[id2], weight, 0.2, 0.2);
            else
                entry.sim = similarity2(sorted_records[id1], sorted_records[id2], record_wt[id1], record_wt[id2], weight, jaccard_threshold, overlap_threshold);
        }
    }

    sort(simpairs.begin(), simpairs.end(), [](const SimPair &p1, const SimPair &p2) {
        return p1.sim > p2.sim;
    });

    /*********** no post-processing *************/
    simpairs.resize(limits);
    for (auto &entry : simpairs)
    {
        int rowid1 = idmap[entry.id1];
        int rowid2 = idmap[entry.id2];

        int pkey1 = record_ids[rowid1];
        int pkey2 = record_ids[rowid2];

        if (pkey1 < pkey2)
            output << pkey1 << "," << pkey2 << endl;
        else
            output << pkey2 << "," << pkey1 << endl;
    }
    /*********** no post-processing *************/

    /***********  post-processing: transitivity  ********/
    /*
    DSU dsu(rec_num);
    for (auto &simpair : simpairs)
    {
        int id1 = simpair.id1;
        int id2 = simpair.id2;
        if (dsu.join(id1, id2))
            if (dsu.over_limit(limits))
                break;
    }
    vector<vector<int>> groups;
    dsu.get_all_groups(groups);

    vector<SimPair> transpairs;
    for (auto &group : groups)
    {
        for (int i = 0; i < group.size(); i++)
        {
            for (int j = i + 1; j < group.size(); j++)
            {
                int id1 = group[i];
                int id2 = group[j];
#ifdef JACCARD_SIMILARITY
                double sim = jaccard(sorted_records[id1], sorted_records[id2]);
#else
                double sim = overlap(sorted_records[id1], sorted_records[id2]);
#endif
                transpairs.emplace_back(id1, id2, sim);
            }
        }
    }

    sort(transpairs.begin(), transpairs.end(), [](const SimPair &p1, const SimPair &p2) {
        return p1.sim > p2.sim;
    });

    transpairs.resize(limits);
    for (auto &entry : transpairs)
    {
        int rowid1 = idmap[entry.id1];
        int rowid2 = idmap[entry.id2];

        int pkey1 = record_ids[rowid1];
        int pkey2 = record_ids[rowid2];

        if (pkey1 < pkey2)
            output << pkey1 << "," << pkey2 << endl;
        else
            output << pkey2 << "," << pkey1 << endl;
    }
    */
    /********** end of post processing *********/
}

void blocking(string &datafilepath, int limits, ofstream &output, 
    double cluster_threshold, 
    bool use_jaccard, double jaccard_threshold, 
    bool use_overlap, double overlap_threshold)
{
    // read the data to tables
    bool normalize = false;
    CSVReader *csvreader = new CSVReader();
    csvreader->reading_one_table(datafilepath, normalize);
    Table &table = csvreader->tables[0];
    // table.PrintInfo();

    /******** get brand info ********/
    vector<string> all_brands;
    all_brands.clear();
    
    int rec_num = table.rows.size();
    vector<int> record_brands(rec_num, -1);

    if (table.schema.size() > 2)
    {
        unordered_map<string, vector<int>> brand_word;
        string delim = " ,-_;.?!|()\\\t\n\r";
        for (int rid = 0; rid < rec_num; rid++)
        {
            string brand = table.rows[rid][3];
            for (int i = 0; i < brand.length(); i++)
                brand[i] = tolower(brand[i]);
            
            vector<string> bag_of_words;
            strToTokens(brand, bag_of_words, delim);
            
            unordered_set<string> bbow;
            for (auto &word : bag_of_words)
                if (word.length() >= 3)
                    bbow.insert(word);

            for (auto &word : bbow)
                brand_word[word].emplace_back(rid);
        }

        for (auto &entry : brand_word)
        {
            // at least 1%
            int minimum_cnt = 0.001 * rec_num;
            if (entry.second.size() > minimum_cnt)
            {
                int brand_id = all_brands.size();
                for (auto tmp : entry.second)
                    record_brands[tmp] = brand_id;
                all_brands.emplace_back(entry.first);
            }
        }
    }
    /******** end brand info ********/

    
    vector<int> record_ids;
    int word_num;
    vector<double> weight;
    vector<int> idmap;
    vector<vector<int>> sorted_records;
    
    txtnorm(table, record_ids, word_num, weight, idmap, sorted_records, all_brands, record_brands);

    // clustering
    vector<vector<int>> cidmap;
    vector<vector<int>> merged_records;

    clustering(cluster_threshold, sorted_records, cidmap, merged_records);

    //similarity join
    postprocessing(record_ids, weight, idmap, cidmap, sorted_records, merged_records, 
        limits, output, use_jaccard, jaccard_threshold, use_overlap, overlap_threshold, record_brands);

    /*
    cout << "first one: " << endl;
    for (auto i = 0; i < sorted_records.front().size(); i++)
        cout << sorted_records.front().at(i) << " ";
    cout << endl;

    cout << "last one: " << endl;
    for (auto i = 0; i < sorted_records.back().size(); i++)
        cout << sorted_records.back().at(i) << " ";
    cout << endl;

    cout << "word num: " << word_num << endl;
    cout << "idmap " << idmap.size() << endl;
    for (auto i = 0; i < idmap.size(); i++)
        cout << idmap[i] << " ";
    cout << endl;
    */
}

int main(int argc, char **argv)
{
    if (argc != 8)
        cout << "datapath, outpath, limit[int], ctheta[0-1], jtheta[0-1], otheta[int], heaprate[int]" << endl;

    string data1 = argv[1];
    string outpath = argv[2];
    int limit = atoi(argv[3]);
    double ctheta = atof(argv[4]);
    double jtheta = atof(argv[5]);
    double otheta = atoi(argv[6]);
    HEAPRATE = atoi(argv[7]);
    
    ofstream out(outpath, ios::out);
    

    blocking(data1, limit, out, ctheta, true, jtheta, true, otheta);  // merging theta, jac theta, overlap theta


    out.close();

    /*
    // calculating recall
    string gt1 = "Y1.csv";
    string gt2 = "Y2.csv";
    string ans = "output.csv";

    ifstream in1(gt1, ios::in);
    ifstream in2(gt2, ios::in);
    ifstream in3(ans, ios::in);

    string str;
    getline(in1, str);
    vector<string> gt1_str;
    while (getline(in1, str))
        gt1_str.emplace_back(str);

    getline(in2, str);
    vector<string> gt2_str;
    while (getline(in2, str))
        gt2_str.emplace_back(str);
    
    getline(in3, str);
    unordered_set<string> ans1_str;
    int lnum = 0;
    while (getline(in3, str))
    {
        ans1_str.insert(str);
        if (++lnum == 100000)
            break;
    }
    int fnum = 0;
    for (auto &entry : gt1_str)
    {
        if (ans1_str.find(entry) != ans1_str.end())
            fnum++;
    }
    cout << fnum << " recall1: " << fnum * 1.0 / gt1_str.size() << endl;

    unordered_set<string> ans2_str;
    while (getline(in3, str))
    {
        ans2_str.insert(str);
    }
    
    fnum = 0;
    for (auto &entry : gt2_str)
    {
        if (ans2_str.find(entry) != ans2_str.end())
            fnum++;
    }
    cout << fnum << " recall2: " << fnum * 1.0 / gt2_str.size() << endl;
    */
}
