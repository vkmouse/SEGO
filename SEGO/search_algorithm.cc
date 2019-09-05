#include "__header__.h"
void cec17_test_func(double *, double *, int, int, int);

void search_algorithm::init(string output_folder, int num_runs, int func_num, int num_dim)
{
    this->output_folder = output_folder;
    this->num_runs = num_runs;
    this->func_num = func_num;
    this->num_dim = num_dim;
    max_nfes = num_dim * 10000;
    optimal = func_num * 100;
    record_obj.resize(num_runs, v1d<double>(num_record, 0.0));
    cout << num_dim << "D\tproblem = " << func_num << endl;
    for (int i = 0; i < num_record; i++)
        record_period[i] = (i + 1) / (double)num_record;
}
void search_algorithm::run_all()
{
    clock_t t1, t2, total_time;
    total_time = 0.0;
    for (curr_run = 0; curr_run < num_runs; curr_run++)
    {
        // set random seed to the current time.
        seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        t1 = clock();

        run();

        t2 = clock();
        total_time += t2 - t1;
        cout << curr_run << '\t' << best_sol.fitness - optimal << '\t' << t2 - t1 << endl;
    }
    record_result(total_time / (double)(CLOCKS_PER_SEC) / num_runs); // record parameters, result, and time.
}
double search_algorithm::evaluate(v1d<double> sol)
{
    double f;
    // void cec17_test_func(double *, double *, int, int, int);
    cec17_test_func(&sol[0], &f, num_dim, 1, func_num); // 1 is only one objective value
    if (f - optimal < 1e-8)
        f = optimal;
    curr_nfes++;

    if(curr_nfes == (int)(max_nfes * record_period[curr_record_index]))
    {
        record_obj[curr_run][curr_record_index] = best_sol.fitness - optimal;
        curr_record_index++;
    }
    
    return f;
}
void search_algorithm::record_result(double exe_time)
{
    double objectvalue = 0.0;
    for (int i = 0; i < num_runs; i++)
        objectvalue += record_obj[i].back();
    objectvalue /= num_runs;
    cout << objectvalue << endl;

    // record file
    fstream file;
    char fn[256];
    sprintf(fn, "%s/%s_%d_%d.txt", output_folder.c_str(), "se", func_num, num_dim);
    file.open(fn, ios::out);
    for (int i = 0; i < num_runs; i++)
        for (int j = 0; j < num_record; j++)
            file << std::to_string((int)(record_period[j] * max_nfes)) << ',' << record_obj[i][j] << endl;
    file.close();
}
bool search_algorithm::update_best_sol(Solution sol)
{
    if (sol.fitness >= best_sol.fitness) // worst than best
        return false;
    best_sol.sol = sol.sol;
    best_sol.fitness = sol.fitness;
    return true;
}