#include <vector>
#include <string>


typedef std::vector< std::vector<float> > mat_t;
typedef std::vector<float> vec_t;
typedef std::vector<int> vec_int_t;

vec_t read_row(std::istream& str);
void read_ucr_file(std::string fname, vec_int_t *labels, mat_t *data);
