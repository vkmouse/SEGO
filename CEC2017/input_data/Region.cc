#include "__header__.h"
se::Region::Region(int num_dim, v1d<double> max_dim, v1d<double> min_dim, 
		int num_sample, v1d<Solution> sample_sol, v1d<double> sample_fitness, unsigned seed)
{
	this->num_dim = num_dim;
	this->max_dim = max_dim;
	this->min_dim = min_dim;

	this->num_sample = num_sample;
	this->sample_sol = sample_sol;
	this->sample_fitness = sample_fitness;

	G.seed(seed);
	U.param(std::uniform_real_distribution<double>::param_type(0.0, 1.0));

	local_best_obj = 10e99;
	for (int k = 0; k < num_sample; k++)
		if (sample_fitness[k] < local_best_obj)
		{
			local_best_sol = sample_sol[k];
			local_best_obj = sample_fitness[k];
		}

	invest_times = 1;
	n_invest_times = 1;
}
void se::Region::set_num_searcher(v1d<Parameter> param)
{
	idx.clear();
	searcher_sol.clear();
	searcher_fitness.clear();
	searcher_param = param;

	num_searcher = searcher_param.size();
	idx.resize(num_searcher);
	searcher_sol.resize(num_searcher);
	searcher_fitness.resize(num_searcher);

	v1d<int> sample_idx(num_sample, 0);
	for (int k = 0; k < num_sample; k++)
		sample_idx[k] = k;
		
	// select sample
	for (int i = 0; i < num_searcher; i++)
	{
		idx[i] = U(G) * sample_idx.size();
		searcher_sol[i] = sample_sol[idx[i]];
		sample_idx.erase(sample_idx.begin() + idx[i]);
		if (sample_idx.size() == 0)
		{
			sample_idx.assign(num_sample, 0);
			for (int k = 0; k < num_sample; k++)
				sample_idx[k] = k;
		}
	}

	// update invest times
	if (num_searcher == 0)
	{
		invest_times = 1;
		n_invest_times++;
	}
	else
	{
		invest_times++;
		n_invest_times = 1;
	}
}
void se::Region::transit()
{
	for (int i = 0; i < num_searcher; i++)
	{
		int strategy = searcher_param[i].st;
		// transit
		if (strategy == 0)
			searcher_sol[i] = current_to_pbest_2(searcher_sol[i], searcher_param[i].CR, searcher_param[i].F);
		else if (strategy == 1)
			searcher_sol[i] = random_to_pbest_2(searcher_sol[i], searcher_param[i].CR, searcher_param[i].F);
	}
}
void se::Region::set_fitness(v1d<double> fitness)
{
    searcher_fitness = fitness;
}
void se::Region::selection()
{
    for (int i = 0; i < num_searcher; i++)
    {
		if (searcher_fitness[i] < sample_fitness[idx[i]])
		{
			sample_sol[idx[i]] = searcher_sol[i];
			sample_fitness[idx[i]] = searcher_fitness[i];

			if (searcher_fitness[i] < local_best_obj)
			{
				local_best_sol = searcher_sol[i];
				local_best_obj = searcher_fitness[i];
			}
		}
    }
}
Solution se::Region::current_to_pbest_2(Solution curr_sol, double CR, double F)
{
	Solution child_sol(num_dim, 0.0);

	int num_player = num_sample * 0.1;

	double *pbest1_sol, *pbest2_sol, *rand_sol;
	int pbest1 = U(G) * num_sample;
	int pbest2 = U(G) * num_sample;
	int random = U(G) * num_sample;

	// swap
	if (sample_fitness[pbest2] < sample_fitness[pbest1])
	{
		int tmp = pbest1;
		pbest1 = pbest2;
		pbest2 = tmp;
	}
	// touranment
	for (int j = 0; j < num_player - 2; j++)
	{
		int r = U(G) * num_sample;
		if (sample_fitness[r] < sample_fitness[pbest1])
			pbest1 = r;
		if (sample_fitness[r] < sample_fitness[pbest2] && r != pbest1)
			pbest2 = r;
	}
	// repeat
    while (random == pbest1 || random == pbest2)
        random = U(G) * num_sample;

	pbest1_sol = &sample_sol[pbest1][0];
	pbest2_sol = &sample_sol[pbest2][0];
	rand_sol = &sample_sol[random][0];

	// curr + (F * pbest1 -curr) + (F * pbest2 - rand)
	int d_rand = U(G) * num_dim;
	for (int d = 0; d < num_dim; d++)
	{
		double f = U(G);
		if (f <= CR || d == d_rand)
		{
			child_sol[d] = curr_sol[d] + 
						   F * (pbest1_sol[d] - curr_sol[d]) +
					       F * (pbest2_sol[d] - rand_sol[d]);

            if (child_sol[d] > max_dim[d])
                child_sol[d] = (max_dim[d] + curr_sol[d]) / 2.0;
            else if (child_sol[d] < min_dim[d])
                child_sol[d] = (min_dim[d] + curr_sol[d]) / 2.0;
		}
		else
			child_sol[d] = curr_sol[d];
	}

	return child_sol;
}
Solution se::Region::random_to_pbest_2(Solution curr_sol, double CR, double F)
{
	Solution child_sol(num_dim, 0.0);

	int num_player = num_sample * 0.1;

	double *pbest1_sol, *pbest2_sol, *rand1_sol, *rand2_sol;
	int pbest1 = U(G) * num_sample;
	int pbest2 = U(G) * num_sample;
	int rand1 = U(G) * num_sample;
    int rand2 = U(G) * num_sample;

	// swap
	if (sample_fitness[pbest2] < sample_fitness[pbest1])
	{
		int tmp = pbest1;
		pbest1 = pbest2;
		pbest2 = tmp;
	}
	// touranment
	for (int j = 0; j < num_player - 2; j++)
	{
		int r = U(G) * num_sample;
		if (sample_fitness[r] < sample_fitness[pbest1])
			pbest1 = r;
		if (sample_fitness[r] < sample_fitness[pbest2] && r != pbest1)
			pbest2 = r;
	}
	// repeat
    while (rand1 == pbest1 || rand1 == pbest2)
        rand1 = U(G) * num_sample;
    while (rand2 == pbest1 || rand2 == pbest2 || rand2 == rand1)
        rand2 = U(G) * num_sample;

	pbest1_sol = &sample_sol[pbest1][0];
	pbest2_sol = &sample_sol[pbest2][0];
    rand1_sol = &sample_sol[rand1][0];
    rand2_sol = &sample_sol[rand2][0];

	// r1 + (F * pbest1 - r1) + (F * pbest2 - r2)
	int d_rand = U(G) * num_dim;
	for (int d = 0; d < num_dim; d++)
	{
		double f = U(G);
		if (f <= CR || d == d_rand)
		{
			child_sol[d] = rand1_sol[d] + 
                           F * (pbest1_sol[d] - rand1_sol[d]) +
                           F * (pbest2_sol[d] - rand2_sol[d]);

            if (child_sol[d] > max_dim[d])
                child_sol[d] = (max_dim[d] + curr_sol[d]) / 2.0;
            else if (child_sol[d] < min_dim[d])
                child_sol[d] = (min_dim[d] + curr_sol[d]) / 2.0;
		}
		else
			child_sol[d] = curr_sol[d];
	}

	return child_sol;
}
void se::Region::get_success_param(v3d<Parameter> &success_param, v3d<double> &success_diff)
{
	for (int i = 0; i < num_searcher; i++)
	{
		if (searcher_fitness[i] < sample_fitness[idx[i]])
		{
			int j = searcher_param[i].st;
			int k = searcher_param[i].idx;
			success_param[j][k].push_back(searcher_param[i]);
			success_diff[j][k].push_back(sample_fitness[idx[i]] - searcher_fitness[i]);
		}
	}
}
double se::Region::get_T()
{
	return n_invest_times / (double)invest_times;
}
double se::Region::get_V()
{
	// mean of the differential fitness
	double sum = 0.0;
	int count = 0;
    for (int i = 0; i < num_searcher; i++)
    {
		if (searcher_fitness[i] < sample_fitness[idx[i]])
		{
			sum += sample_fitness[idx[i]] - searcher_fitness[i];
			count++;
		}
	}
	if (count == 0)
		return 0;
	return sum / count;
}
double se::Region::get_M()
{
	// local best fitness
	return local_best_obj;
}
void se::Region::reset_num_sample(int num_sample)
{
	while (this->num_sample > num_sample)
	{
		int remove_index = 0;
		for (int k = 1; k < this->num_sample; k++)
			if (sample_fitness[remove_index] < sample_fitness[k])
				remove_index = k;

		sample_sol.erase(sample_sol.begin() + remove_index);
		sample_fitness.erase(sample_fitness.begin() + remove_index);
		this->num_sample--;
	}
}
void se::Region::set_max_min_dim(v1d<double> max_dim, v1d<double> min_dim, v1d<Solution> sol, v1d<double> fitness)
{
	this->max_dim = max_dim;
	this->min_dim = min_dim;
	invest_times = 1;
	n_invest_times = 1;
	sample_sol = sol;
	sample_fitness = fitness;
}
void se::Region::get_sample(v1d<Solution> &sol, v1d<double> &fitness)
{
	sol.insert(sol.end(), sample_sol.begin(), sample_sol.end());
	fitness.insert(fitness.end(), sample_fitness.begin(), sample_fitness.end());
}