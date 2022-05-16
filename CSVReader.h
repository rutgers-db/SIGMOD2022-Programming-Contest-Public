
#ifndef _CSV_READER_H_
#define _CSV_READER_H_

#include <istream>
#include <sstream>
#include <cstring>
#include <dirent.h>
#include <sys/stat.h>
#include <string.h>
#include <stdio.h>

#include "Table.h"

class CSVReader {
private:
  char field_delimiter = ','; // '\t'
  int filesize(const char* filename);
  bool ends_with(std::string const & value, std::string const & ending);
  void csv_read_row(std::istream &in, std::vector<std::string> &row, bool isNorm = true);
  bool get_table(const string &filepath, vector<string> &headers, vector<unordered_map<string, int>> &columns, vector<vector<string>> &rows, bool normalize);
  int max_val_len;
public:
  CSVReader() {}
  vector<Table> tables;

  void strNormalize(string &s); // also for the use of query normalization
  bool reading(string &datafilespath, bool normalize);
  void write_one_table(const Table &table, const string &outfilename);
  bool reading_one_table(string &datafilepath, bool normalize);
  int get_max_val_len() { return max_val_len; };
};

#endif
