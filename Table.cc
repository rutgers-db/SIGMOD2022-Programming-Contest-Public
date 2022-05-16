
#include "Table.h"

Table::Table(int id, const string& name) {
  tid = id;
  table_name = name;
}

Table::Table(int id, const string &name, vector<string> &data_headers, vector<vector<string>> &data_rows, vector<unordered_map<string, int>> &data_columns) {
  tid = id;
  table_name = name;

  schema = data_headers;
  rows = data_rows;
  cols = data_columns;

  row_no = rows.size();
  col_no = cols.size();
}

void Table::Profile() {
  row_no = rows.size();
  col_no = cols.size();
}

void Table::PrintInfo() {
  cout << " number of rows: " << row_no << endl;
  cout << " number of columns: " << col_no << endl;
  cout << " the schema: ";
  for (auto &attr : schema) cout << attr << "; ";
  cout << endl;
}