#ifndef _HEADER_H_
#define _HEADER_H_

#include <vector>
#include <chrono>
#include <random>
#include <fstream>
#include <iostream>
#include <algorithm>

using std::cout;
using std::endl;
using std::fstream;
using std::ios;
using std::string;
using std::vector;

struct Solution
{
    vector<double> sol;
    double fitness;
};
template <typename T>
using v1d = vector<T>;
template <typename T>
using v2d = vector<vector<T>>;
template <typename T>
using v3d = vector<vector<vector<T>>>;

class search_algorithm
{
public:
    virtual void run() = 0;
    void init(string output_folder, int num_runs, int func_num, int num_dim);
    void run_all();

    int num_runs;
    int func_num;
    int num_dim;
    int curr_nfes;  // modify: init(),
    int max_nfes;   // maximal number of function evaluations, modify: curr_nfes
    double optimal;
    double max_dim = 100.0;
    double min_dim = -100.0;
    int curr_run;
    Solution best_sol;// modify: init(),
    int num_record = 100;
    int curr_record_index;// modify: init(), curr_record_index()
    double record_period[100];
    v2d<double> record_obj;
    string output_folder;

    // parameter for random instruction
    unsigned seed;
    std::default_random_engine generator;   // modify: init(),
    std::uniform_real_distribution<double> uniform; // modify: init(),

    // method
    double evaluate(v1d<double> sol);
    void record_result(double exe_time);
    bool update_best_sol(Solution sol);
};
class se : public search_algorithm
{
public:
    virtual void run();
    se(double init_searcher, int init_region, double sample_rate, int memory_size, double touranment_rate, double rd2);

private:
    // input parameter, only assign in the constuction.
    double searcher_rate;   // alpha
    int init_region;        // beta
    double sample_rate;     // gamma
    int memory_size;        // H
    double touranment_rate; // t
    double rd2;             // rd

    typedef struct Parameter
    {
        double CR;
        double F;
        int strategy;       // ST
    } Parameter;
    typedef struct Searcher
    {
        Parameter param;    // modify: assign_parameter()
        Solution trial_vec; // modify: transition(), compute_expected_value()
        int region_index;   // modify: transition()
        int sample_index;   // modify: transition()
    } Searcher;
    typedef struct Region
    {
        v1d<double> max_dim;    // [num_dim], modify: init(), resource_arrangement(), combine_space()
        v1d<double> min_dim;    // [num_dim], modify: init(), resource_arrangement(), combine_space()
        int num_searcher;       // modify: init(), resource_arrangement(), determination(), combine_space()
        int invest_times;       // modify: init(), resource_arrangement(), market_research(), combine_space()
        int not_invest_times;   // modify: init(), market_research(), combine_space()
    } Region;

    int init_searcher;  // modify: init()
    int curr_searcher;  // modify: init(), linear_reduction()
    int num_region;     // modify: init(), linear_reduction(), combine_space()
    int num_sample;     // modify: init(), linear_reduction(), combine_space()
    int min_searcher = 4;
    int num_strategy = 2;

    v1d<Region> regions;        // [num_region], modify: init(), determination(), market_research(), combine_space()
    v2d<Solution> samples;      // [num_region][num_sample], modify: init(), resource_arrangement(), update_parameter(), linear_reduction(), combine_space()
    v1d<Searcher> searchers;    // [num_searcher], modify: init(), assign_parameter(), transition(), compute_expected_value(), linear_reduction()
    v2d<Parameter> memory_param;// [num_strategy][memory_size], modify: init(), update_parameter(),
    int memory_position;        // [num_strategy], modify: init(), update_parameter(),
    v1d<double> expected_value; // [num_region], modify: init(), combine_space()
    v1d<int> split_dim;         // [log(num_region)], modify: resource_arrangement(), combine_space()
    std::normal_distribution<double> normal;    // modify: assign_parameter()
    std::cauchy_distribution<double> cauchy;    // modify: assign_parameter()

    void init();
    v1d<Searcher> transition(v1d<Region>, v2d<Solution>);
    v1d<Searcher> assign_parameter(v2d<Parameter>);
    v1d<double> current_to_tbest(int, int, Parameter, int);
    v1d<double> random_to_tbest(int, Parameter, int);
    v2d<Parameter> update_parameter(v1d<Searcher>);
    void devide_space();
    v1d<double> compute_expected_value(v2d<Solution>);
    v1d<Region> determination(v1d<double>);
    v1d<double> normalization(v1d<double>);
    void linear_reduction();
    void destroy();
    void combine_space();
    void resource_arrangement();
    v1d<Region> market_research(v1d<Region>);
};
#endif
