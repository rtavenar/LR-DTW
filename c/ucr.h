#include <vector>
#include <string>

#include "types.h"

vec_t read_row(std::istream& str);
void read_ucr_file(std::string fname, vec_int_t *labels, mat3d_t *data);
