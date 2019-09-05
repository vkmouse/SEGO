#include "__header__.h"

double *OShift, *M, *y, *z, *x_bound;
int ini_flag = 0, n_flag, func_flag, *SS;

int main(int argc, char **argv)
{
    search_algorithm *search_alg;

    search_alg = new se(
        atof(argv[5]), //max_searcher
        atoi(argv[6]), //max_region
        atof(argv[7]), //sample_rate
        atoi(argv[8]), //memory_size
        atof(argv[9]), //touranment_rate
        atof(argv[10]) //rd
    );

    search_alg->init(argv[1],
                     atoi(argv[2]),
                     atoi(argv[3]),
                     atoi(argv[4]));
    search_alg->run_all();
    return 0;
}