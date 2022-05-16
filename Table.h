
#ifndef _TABLE_H_
#define _TABLE_H_

#include <iostream>
#include <fstream>
#include <vector> 
#include <string>
#include <unordered_map>

using namespace std;

class Table {
public:
  int tid;
  int row_no, col_no;
  string table_name;
  vector<string> schema;
  vector<vector<string>> rows;
  vector<unordered_map<string, int>> cols;
  
  Table (int id, const string &name);
  Table (int id, const string &name, vector<string> &data_headers, vector<vector<string>> &data_rows, 
      vector<unordered_map<string, int>> &data_columns);
  void Profile();
  void PrintInfo();
};
#endif
