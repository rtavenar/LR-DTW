#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "ucr.h"



vec_t read_row(std::istream& str) {
    vec_t result;
    std::string line, cell;

    std::getline(str, line);
    std::stringstream line_stream(line);

    while(std::getline(line_stream, cell, ' ')) {
        if(cell.size() > 0)
            result.push_back(std::stof(cell));
    }
    return result;
}

void read_ucr_file(std::string fname, vec_int_t *labels, mat_t *data) {
    std::ifstream file;
    file.open(fname);
    vec_t row;

    row = read_row(file);
    while (row.size() > 0) {
        (*labels).push_back(static_cast<int>(row[0]));
        vec_t row_data;
        for(size_t i=1; i<row.size(); ++i)
            row_data.push_back(row[i]);
        (*data).push_back(row_data);
        row = read_row(file);
    }
}
