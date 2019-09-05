#include "__header__.h"
void cec17_test_func(double *x, double *f, int nx, int mx,int func_num);
se::se(double searcher_rate,int init_region, double sample_rate, int memory_size, double touranment_rate, double rd2)
{
	this->searcher_rate = searcher_rate;
	this->init_region = init_region;
	this->sample_rate = sample_rate;
	this->memory_size = memory_size;
	this->touranment_rate = touranment_rate;
	this->rd2 = rd2;
	
	cout << "searcher_rate = " << searcher_rate << endl;
	cout << "init_region = " << init_region << endl;
	cout << "sample_rate = " << sample_rate << endl;
	cout << "memory_size = " << memory_size << endl;
	cout << "touranment_rate = " << touranment_rate << endl;
	cout << "rd2 = " << rd2 << endl;
}
void se::run()
{
	init();
	resource_arrangement();
	while (curr_nfes < max_nfes)
	{
		(v1d<Searcher>)searchers = assign_parameter(memory_param);
		(v1d<Searcher>)searchers = transition(regions, samples);// ----	|
		expected_value = compute_expected_value(samples);		//		|--	vision search
		regions = determination(expected_value);				// ----	|
		regions = market_research(regions);
		(v2d<Parameter>)memory_param = update_parameter(searchers);
		linear_reduction();
		for (int i = 0; i < curr_searcher; i++)
			update_best_sol(searchers[i].trial_vec);
	}
	for (int j = 0; j < num_dim; j++)
	{
		if (best_sol.sol[j] > max_dim || best_sol.sol[j] < min_dim)
		{
			for (int i = 0; i < num_dim - 1; i++)
				cout << best_sol.sol[i] << ", ";
			cout << best_sol.sol.back() << endl;
			break;
		}
	}
	destroy();
}

/** @init(): Initialize the variables
 * 	1. 	Initialize the variables which is used by search algorithm.
 * 	2. 	Initialize the variables which is used by SE and assign the memory size for the vector.
 * 	3. 	Initialize the parameter's memory including CR = 0.5, F = 0.5.
 *  	There are two parameter's memories. 
 * 	  	One is for current-to-tbest, the other is for random-to-tbest.
 * 
 * 	param in: 	null
 *	param out: 	null
 * 	unmodified golbal value:
 * 				searcher_rate, init_region, sample_rate, 
 * 				num_dim, max_dim, min_dim, 
 * 				num_strategy, num_memory
 * 	modified golbal value:
 * 				curr_nfes, generator, uniform, best_sol, curr_record_index (search algorithm variables)
 * 				init_searcher, curr_searcher, num_region, num_sample, searchers, regions, samples, 
 * 				memory_param, memory_position, expected_value
*/
void se::init()
{
	// 1. initialize search algorithm variables
	curr_nfes = 0;
	generator.seed(seed);
	uniform.param(std::uniform_real_distribution<double>::param_type(0.0, 1.0));
	best_sol = Solution{.sol = v1d<double>(num_dim, 0.0), .fitness = 10e99};
	curr_record_index = 0;

	// 2. initalize SE variables
	init_searcher = num_dim * searcher_rate;
	curr_searcher = init_searcher;
	num_region = init_region;
	num_sample = (int)round(curr_searcher / (double)num_region * sample_rate);

	searchers.resize(curr_searcher);
	regions.assign(num_region, Region{
					.max_dim = v1d<double>(num_dim, max_dim), 
					.min_dim = v1d<double>(num_dim, min_dim),
					.num_searcher = curr_searcher / num_region,
					.invest_times = 1,
					.not_invest_times = 1});
	samples.assign(num_region, v1d<Solution>(num_sample, Solution{.sol = v1d<double>(num_dim, 0.0), .fitness = 0.0}));
	memory_param.assign(num_strategy, v1d<Parameter>(memory_size, Parameter{}));
	memory_position = 0;
	expected_value.resize(num_region);

	// 3. initialize memory which is used to adaptive CR and F.
	for (int k = 0; k < num_strategy; k++)
		for (int h = 0; h < memory_size; h++)
		{
			memory_param[k][h].CR = 0.5;
			memory_param[k][h].F = 0.5;
			memory_param[k][h].strategy = k;
		}
}
/** @resource_arrangement(): Divide the search space and assign computation resources to each region.
 *	1.	Divide the search space and set the limitation boundary to each region.
 *	2.	Generate the population to each region.
 *	3.	Equaly assign computation resources (number of searchers) to each region.
 *
 * 	param in: 	null
 *	param out: 	null
 * 	unmodified golbal value:
 * 				num_region, num_dim, max_dim, min_dim, num_sample
 * 	modified golbal value:
 * 				split_dim, regions, samples
*/
void se::resource_arrangement()
{
	// 1. devide search space
	// number of dimension should be splited
	int num_split = log2(num_region); 

	// randomly select num_split numbers to be divided dimensions and store them.
	for (int i = 0; i < num_split; i++)
	{
		int r;
		do
		{
			r = uniform(generator) * num_dim;
		} while (std::find(split_dim.begin(), split_dim.end(), r) != split_dim.end());
		split_dim.push_back(r);
	}

	// convert region index from decimal to binary.
	v2d<bool> bin_region_index(num_region, v1d<bool>(num_split, 0));
	for (int j = 0; j < num_region; j++)
	{
		int region_index = j;
		for (int i = 0; i < num_split; i++)
		{
			if (region_index % 2) // remainder is 1
				bin_region_index[j][i] = 1;
			region_index /= 2;
		}
	}

	// resize the boundary of each region
	for (int j = 0; j < num_region; j++)
	{
		for (int i = 0; i < num_split; i++)
		{
			if (bin_region_index[j][i]) // is one
				// set the boundary of dimension split_dim[i] to front half.
				// but also keep 1/10 overlap.
				regions[j].max_dim[split_dim[i]] = (max_dim + min_dim) / 2.0 + (max_dim - min_dim) / 10.0; // se
			else
				regions[j].min_dim[split_dim[i]] = (max_dim + min_dim) / 2.0 - (max_dim - min_dim) / 10.0; // se
		}
	}

	// 2. initialize region's population
	for (int j = 0; j < num_region; j++)
		for (int k = 0; k < num_sample; k++)
		{
			for (int d = 0; d < num_dim; d++)
				samples[j][k].sol[d] = uniform(generator) * (regions[j].max_dim[d] - regions[j].min_dim[d]) + regions[j].min_dim[d];
			samples[j][k].fitness = evaluate(samples[j][k].sol);
			update_best_sol(samples[j][k]);
		}

	// 3. initialize computation
	int sum = 0;
	for (Region r: regions)
		sum += r.num_searcher;
	for (int j = 0; sum != curr_searcher; j = (j + 1) % num_region)
	{
		regions[j].num_searcher++;
		sum++;
		regions[j].invest_times = regions[j].num_searcher;
	}
}
/** @assign_parameter(): Assign the search parameter to each searchers.
 * 	1.	Select a following parameter.
 * 	2.	Generate a new CR.
 * 	3.	Generate a new F.
 * 
 * 	param in: 	memory_param
 *	param out: 	searchers
 * 	unmodified golbal value:
 * 				memory_size, num_strategy
 * 	modified golbal value:
 * 				null
*/
v1d<se::Searcher> se::assign_parameter(v2d<Parameter> memory_param) // Randomly select from 1 to memory size and record the index
{
	for (Searcher &searcher : searchers)
	{
		// 1. Randomly select a transition strategy and a parameter in the memory.
		int index = uniform(generator) * memory_size;
		int strategy = searcher.param.strategy = uniform(generator) * num_strategy;

		// 2. generate a CR and repairs its value.
		// memory's CR is -1 when improving failed in the previous iteration.
		if (memory_param[strategy][index].CR == -1)
			searcher.param.CR = 0;
		else
		{
			normal.param(std::normal_distribution<double>::param_type(memory_param[strategy][index].CR, 0.1));
			searcher.param.CR = normal(generator);
			if (searcher.param.CR > 1.0)
				searcher.param.CR = 1.0;
			else if (searcher.param.CR < 0.0)
				searcher.param.CR = 0.0;
		}

		// 3. generate a F and repair its value
		// only a little change in the early period.
		if (max_nfes * rd2 > curr_nfes)
		{
			cauchy.param(std::cauchy_distribution<double>::param_type(memory_param[strategy][index].F, 0.01));
			searcher.param.F = cauchy(generator);
			if (searcher.param.F > 0.6)
				searcher.param.F = 0.6;
			else if (searcher.param.F < 0.5)
				searcher.param.F = 0.5;
		}
		else
		{
			do 
			{
				cauchy.param(std::cauchy_distribution<double>::param_type(memory_param[strategy][index].F, 0.1));
				searcher.param.F = cauchy(generator);
			} while (searcher.param.F <= 0);
			if (searcher.param.F > 1)
				searcher.param.F = 1;
		}
	}
	return searchers;
}
/**	@transition(): Move the solution according the strategy.
 * 	1.	Assign a sample to each searchers.
 * 	2.	Calculate the tournament times.
 * 	3.	Execute the transition strategy (current-to-tbest or random-to-tbest)
 * 
 * 	param in: 	regions, samples
 *	param out: 	searchers
 * 	unmodified golbal value:
 * 				num_region, num_sample, touranment_rate, curr_searcher
 * 	modified golbal value:
 * 				null
*/
v1d<se::Searcher> se::transition(v1d<Region> regions, v2d<Solution> samples)
{
	// 1. Select a sample solution being the comparison solution
	// Each region will execute the transition in its sample population.
	for (int i = 0, j = 0; j < num_region; j++)
	{
		// the number of computation resource invested in this region.
		int num_searcher = regions[j].num_searcher;

		// Every sample will be assigned to the searcher,
		// while the number of searchers is greater than the number of samples.
		while (num_searcher >= num_sample)
		{
			for (int m = 0; m < num_sample; m++)
			{
				searchers[i].region_index = j;
				searchers[i].sample_index = m;
				i++;
			}
			num_searcher -= num_sample;
		}
		// if 0 < num_searcher < num_sample, randomly select unique samples.
		if (num_searcher != 0)
		{
			v1d<int> index(num_sample);
			for (int k = 0; k < num_sample; k++)
				index[k] = k;
			shuffle (index.begin(), index.end(), generator);
			for (int m = 0; m < num_searcher; m++)
			{
				searchers[i].region_index = j;
				searchers[i].sample_index = index[m];
				i++;
			}
		}
	}

	// 2. set tournament times
	// tournament times set to num_sample multiply by touranment_rate
	int num_player;
	// tournament times will be double in the final period
	if (curr_nfes < 0.8 * max_nfes)
	 	num_player = round(num_sample * touranment_rate);	
	else
		num_player = round(num_sample * touranment_rate * 2.0);
	if (num_player < 2)
		num_player = 2;
	
	// 3. execute the transit
	for (int i = 0; i < curr_searcher; i++)
	{
		if (searchers[i].param.strategy == 0)
			searchers[i].trial_vec.sol = current_to_tbest(searchers[i].region_index, searchers[i].sample_index, searchers[i].param, num_player);
		else if (searchers[i].param.strategy == 1)
			searchers[i].trial_vec.sol = random_to_tbest(searchers[i].region_index, searchers[i].param, num_player);
	}
	return searchers;
}
/** @current_to_tbest(): Move the vector according current solution and tbest.
 * 	1.	Select the following tbest solution.
 * 	2.	Move the vector in the boundary.
 * 	
 * 	param in: 	region_index, current_solution_index, param, num_player
 *	param out: 	child_sol (only a solution)
 * 	unmodified golbal value:
 * 				num_sample, samples, num_dim, regions,
 * 	modified golbal value:
 * 				null
*/
v1d<double> se::current_to_tbest(int region_index, int curr, Parameter param, int num_player)
{
	// 1. Select the following solution
	int tbest1 = uniform(generator) * num_sample;
	int tbest2 = uniform(generator) * num_sample;
	int random = uniform(generator) * num_sample;

	// generate unique following solution
    while (tbest1 == tbest2)
		tbest2 = uniform(generator) * num_sample;
	// touranment selection
	for (int j = 0; j < num_player - 1; j++)
	{
		int r = uniform(generator) * num_sample;
		if (samples[region_index][r].fitness < samples[region_index][tbest1].fitness && r != curr)
			tbest1 = r;
		if (samples[region_index][r].fitness < samples[region_index][tbest2].fitness && r != tbest1 && r != curr)
			tbest2 = r;
	}
    while (random == tbest2 || random == curr)
        random = uniform(generator) * num_sample;

	// convert the solution index to the solution address
	double *curr_sol = &samples[region_index][curr].sol[0];
	double *tbest1_sol = &samples[region_index][tbest1].sol[0];
	double *tbest2_sol = &samples[region_index][tbest2].sol[0];
	double *rand_sol = &samples[region_index][random].sol[0];

	// 2. Move the vector and keep in the limitation boundary.
	v1d<double> child_sol(num_dim, 0.0);
	int d_rand = uniform(generator) * num_dim;
	for (int d = 0; d < num_dim; d++)
	{
		double f = uniform(generator);
		if (f <= param.CR || d == d_rand)
		{
			// curr + (F * tbest1 - curr) + (F * tbest2 - rand)
			child_sol[d] = curr_sol[d] + 
						   param.F * (tbest1_sol[d] - curr_sol[d]) +
					       param.F * (tbest2_sol[d] - rand_sol[d]);

			// keep in the boundary
            if (child_sol[d] > regions[region_index].max_dim[d])
                child_sol[d] = (regions[region_index].max_dim[d] + curr_sol[d]) / 2.0;
            else if (child_sol[d] < regions[region_index].min_dim[d])
                child_sol[d] = (regions[region_index].min_dim[d] + curr_sol[d]) / 2.0;
		}
		else
			child_sol[d] = curr_sol[d];
	}

	return child_sol;
}
/** @random_to_tbest(): Move the vector according random solution and tbest.
 * 	1.	Select the following tbest solution.
 * 	2.	Move the vector in the boundary.
 * 	
 * 	param in: 	region_index, param, num_player
 *	param out: 	child_sol (only a solution)
 * 	unmodified golbal value:
 * 				num_sample, samples, num_dim, regions,
 * 	modified golbal value:
 * 				null
*/
v1d<double> se::random_to_tbest(int region_index, Parameter param, int num_player)
{
	// 1. Select the following solution
	int tbest1 = uniform(generator) * num_sample;
	int tbest2 = uniform(generator) * num_sample;
	int rand1 = uniform(generator) * num_sample;
    int rand2 = uniform(generator) * num_sample;

	// generate unique following solution
    while (tbest1 == tbest2)
		tbest2 = uniform(generator) * num_sample;
	// touranment selection
	for (int j = 0; j < num_player - 1; j++)
	{
		int r = uniform(generator) * num_sample;
		if (samples[region_index][r].fitness < samples[region_index][tbest1].fitness && r != rand1)
			tbest1 = r;
		if (samples[region_index][r].fitness < samples[region_index][tbest2].fitness && r != tbest1 && r != rand1)
			tbest2 = r;
	}
    while (rand2 == tbest2 || rand2 == rand1)
        rand2 = uniform(generator) * num_sample;

	// convert the solution index to the solution address
	double *tbest1_sol = &samples[region_index][tbest1].sol[0];
	double *tbest2_sol = &samples[region_index][tbest2].sol[0];
	double *rand1_sol = &samples[region_index][rand1].sol[0];
	double *rand2_sol = &samples[region_index][rand2].sol[0];

	// 2. Move the vector and keep in the limitation boundary.
	v1d<double> child_sol(num_dim, 0.0);
	int d_rand = uniform(generator) * num_dim;
	for (int d = 0; d < num_dim; d++)
	{
		double f = uniform(generator);
		if (f <= param.CR || d == d_rand)
		{
			// rand1 + (F * tbest1 - rand1) + (F * tbest2 - rand2)
			child_sol[d] = rand1_sol[d] + 
						   param.F * (tbest1_sol[d] - rand1_sol[d]) +
					       param.F * (tbest2_sol[d] - rand2_sol[d]);

			// keep in the boundary
            if (child_sol[d] > regions[region_index].max_dim[d])
                child_sol[d] = (regions[region_index].max_dim[d] + rand1_sol[d]) / 2.0;
            else if (child_sol[d] < regions[region_index].min_dim[d])
                child_sol[d] = (regions[region_index].min_dim[d] + rand1_sol[d]) / 2.0;
		}
		else
			child_sol[d] = rand1_sol[d];
	}

	return child_sol;
}
/**	@compute_expected_value(): Compute expected value depend on three factors.
 * 	1.	evaluate the searchers which does not have fitness after transition operator.
 * 	2.	compute three influence factors of expected value.
 * 	3.	normalize the factors and compute the expected value.
 * 
 * 	param in: 	samples
 *	param out: 	expected_value
 * 	unmodified golbal value:
 * 				curr_searcher, num_region, regions
 * 	modified golbal value:
 * 				searchers
*/
v1d<double> se::compute_expected_value(v2d<Solution> samples)
{
	// 1. evaluation
	for (int i = 0; i < curr_searcher; i++)
		searchers[i].trial_vec.fitness = evaluate(searchers[i].trial_vec.sol);

	// 2. Compute three factor T, V, and M respectively.
	v1d<double> T(num_region), V(num_region), M(num_region);
	v1d<double> sum(num_region, 0.0);
	v1d<int> count(num_region, 0);
	// sum the improving fitness.
	for (int i = 0; i < curr_searcher; i++)
	{
		int j = searchers[i].region_index;
		int k = searchers[i].sample_index;
		if (searchers[i].trial_vec.fitness < samples[j][k].fitness)
		{
			sum[j] += samples[j][k].fitness - searchers[i].trial_vec.fitness;
			count[j]++;
		}
	}
	// calculate the factors of each region.
	for (int j = 0; j < num_region; j++)
	{
		// T is investment situation
		T[j] = regions[j].not_invest_times/ (double)regions[j].invest_times;
		// V is average improvement
		if (count[j] == 0)
			V[j] = 0;
		else
			V[j] = sum[j] / count[j]; 
		// M is the best fitness in the region.
		M[j] = samples[j][0].fitness;
		for (int k = 1; k < num_sample; k++)
			if (samples[j][k].fitness < M[j])
				M[j] = samples[j][k].fitness;
	}

	// 3. normalize and calculate expected value
	T = normalization(T);
	V = normalization(V);
	M = normalization(M);
	for (int j = 0; j < num_region; j++)
	{
		T[j] += 0.001;
		V[j] += 0.001;
		M[j] = 1.0 - M[j] + 0.001; // greater is better
	}
	// multiply three factors
	for (int j = 0; j < num_region; j++)
		expected_value[j] = T[j] * V[j] * M[j];
	
	return expected_value;
}
/**	@determination(): Determine the number of computation resource should be invested in each region.
 * 	1.	Initialize the number of resource computation.
 * 	2.	Searcher chooses the region using tournament selection.
 * 
 * 	param in: 	expected_value
 *	param out: 	regions
 * 	unmodified golbal value:
 * 				curr_searcher, num_region
 * 	modified golbal value:
 * 				null
*/
v1d<se::Region> se::determination(v1d<double> expected_value)
{
	// 1. Initialize the number of searcher of each region.
	for (int j = 0; j < num_region; j++)
		regions[j].num_searcher = 0;

	// 2. Select the region of each searcher in next generation.
	// tournament times set to half number of region.
	int num_player = num_region / 2;
	for (int i = 0; i < curr_searcher; i++)
	{
		int selected_region = uniform(generator) * num_region;
		for (int k = 1; k < num_player; k++)
		{
			int rand_region = uniform(generator) * num_region;
			if (expected_value[selected_region] < expected_value[rand_region])
				selected_region = rand_region;
		}
		regions[selected_region].num_searcher++;
	}

	return regions;
}
/**	@market_research(): Update the invested information of each region.
 * 	1.	Update invested information
 * 
 * 	param in: 	regions
 *	param out: 	regions
 * 	unmodified golbal value:
 * 				num_region
 * 	modified golbal value:
 * 				null
*/
v1d<se::Region> se::market_research(v1d<Region> regions)
{
	// 1. update invest status
	for (int j = 0; j < num_region; j++)
	{
		int num_searcher = regions[j].num_searcher;

		if (num_searcher == 0)
		{
			regions[j].invest_times = 1;
			regions[j].not_invest_times++;
		}
		else
		{
			regions[j].invest_times += num_searcher;
			regions[j].not_invest_times = 1;
		}
	}
	return regions;
}
/**	@update_parameter(): Update the memory parameter depend on the parameter improving successfully.
 * 	1.	Collect the successful parameter and its improvement.
 * 	2.	Update the sample solution and fitness based on the comparison fitness of searcher.
 * 	3.	Update the pos-th memory with weight lehmer mean.
 * 
 * 	param in: 	searchers
 *	param out: 	memory_param
 * 	unmodified golbal value:
 * 				curr_searcher, num_strategy
 * 	modified golbal value:
 * 				samples, memory_position
*/
v2d<se::Parameter> se::update_parameter(v1d<Searcher> searchers)
{
	// Record the parameter improving successfully.
	v2d<Parameter> success_param(num_strategy);
	// Record the improvement of successful parameter
	v2d<double> success_diff(num_strategy);

	// 1. Collect the information of successful improvement.
	for (int i = 0; i < curr_searcher; i++)
	{
		int j = searchers[i].region_index;
		int k = searchers[i].sample_index;
		if (searchers[i].trial_vec.fitness < samples[j][k].fitness)
		// children is better than origin
		// record the successful parameter
		{
			int strategy = searchers[i].param.strategy;
			success_param[strategy].push_back(searchers[i].param);
			success_diff[strategy].push_back(samples[j][k].fitness - searchers[i].trial_vec.fitness);
		}
	}

	// 2. Update sample
	for (int i = 0; i < curr_searcher; i++)
	{
		int j = searchers[i].region_index;
		int k = searchers[i].sample_index;
		if (searchers[i].trial_vec.fitness == samples[j][k].fitness)
			samples[j][k].sol = searchers[i].trial_vec.sol;
		else if (searchers[i].trial_vec.fitness < samples[j][k].fitness)
		// children is better than origin
		{
			samples[j][k].sol = searchers[i].trial_vec.sol;
			samples[j][k].fitness = searchers[i].trial_vec.fitness;
		}
	}

	// 3. Update memory parameter
	for (int k = 0; k < num_strategy; k++)
	{
		int num_success_params = success_param[k].size();
		if (num_success_params > 0)
		{
			double sum_F = 0.0;
			double sum_CR = 0.0;
			double sum_diff = 0.0;
			double sum_pow_CR = 0.0;
			double sum_pow_F = 0.0;

			for (int l = 0; l < num_success_params; l++)
				sum_diff += success_diff[k][l];

			// weight lehmer mean
			for (int l = 0; l < num_success_params; l++)
			{
				double weight = success_diff[k][l] / sum_diff;

				sum_F += weight * success_param[k][l].F;
				sum_pow_F += weight * success_param[k][l].F * success_param[k][l].F;
				sum_CR += weight * success_param[k][l].CR;
				sum_pow_CR += weight * success_param[k][l].CR * success_param[k][l].CR;
			}
			if (sum_F == 0.0)
				memory_param[k][memory_position].F = 0.0;
			else
				memory_param[k][memory_position].F = sum_pow_F / sum_F;
			if (sum_CR == 0.0 || memory_param[k][memory_position].CR == -1)
				memory_param[k][memory_position].CR == -1;
			else
				memory_param[k][memory_position].CR = sum_pow_CR / sum_CR;
		}
	}
	memory_position++;
	if (memory_position >= memory_size)
		memory_position = 0;

	return memory_param;
}
/**	@normalization(): Normalize the array to the range [0,1].
 * 	1.	Find the maxium and minium value.
 * 	2.	Normalize
 * 
 * 	param in: 	array
 *	param out: 	array
 * 	unmodified golbal value:
 * 				null
 * 	modified golbal value:
 * 				null
*/
v1d<double> se::normalization(v1d<double> array)
{
	// 1. find the maxium and minium value.
	double min = array[0], max = array[0];
	for (int i = 1; i < array.size(); i++)
	{
		if (min > array[i])
			min = array[i];
		if (array[i] > max)
			max = array[i];
	}
	// 2. Normalize
	if (max - min == 0.0)
		for (int i = 0; i < array.size(); i++)
			array[i] = 1 / array.size();
	else
		for (int i = 0; i < array.size(); i++)
			array[i] = (array[i] - min) / (max - min);
	return array;
}
/** @linear_reduction(): Reduce the searcher, samples, region.
 * 	1.	Reduce the region.
 * 	2.	Reduce the searchers and samples
 * 
 * 	param in: 	null
 *	param out: 	null
 * 	unmodified golbal value:
 * 				init_region, curr_nfes, max_nfes, rd2, init_searcher,
 * 				min_searcher,
 * 	modified golbal value:
 * 				curr_searcher, searchers, num_sample, samples, num_region
*/
void se::linear_reduction()
{
	int plan_region = init_region / pow(2, (int)(curr_nfes / (max_nfes * rd2 / log2(init_region))));
	if (plan_region < 1)
		plan_region = 1;
	if (num_region > plan_region)
	{
		combine_space();
	}

	int plan_curr_searcher = round(((min_searcher - init_searcher) / (double)max_nfes * curr_nfes) + init_searcher);
	if (curr_searcher > plan_curr_searcher)
	{
		for (int i = 0; i < curr_searcher - plan_curr_searcher; i++)
			regions[uniform(generator) * num_region].num_searcher--;
		curr_searcher = plan_curr_searcher;
		if (curr_searcher < min_searcher)
			curr_searcher = min_searcher;
		searchers.resize(curr_searcher);		// it will change when linear reudution.

		int reduction_num = num_sample - (int)round(curr_searcher / (double)num_region * sample_rate);
		if (num_sample - reduction_num < min_searcher)
			reduction_num = num_sample - min_searcher;

		// reduce population with sort
		for (int i = 0; i < reduction_num; i++)
		{
			for (int j = 0; j < num_region; j++)
			{
				int index = 0;
				for (int k = 0; k < num_sample; k++)
					if (samples[j][k].fitness > samples[j][index].fitness) // worst than index
						index = k;
				samples[j].erase(samples[j].begin() + index);
			}
			num_sample--;
		}
	}
}
/**	@combine_space(): Combine the number of search spaces from n to n/2.
 * 	1.	Combine the region.
 * 	2.	Clear the variables.
 * 
 * 	param in: 	null
 *	param out: 	null
 * 	unmodified golbal value:
 * 				null
 * 	modified golbal value:
 * 				regions, samples, 
 * 				expected_value, num_region, split_dim
*/
void se::combine_space()
{
	// 1. Combine the region.
	int num_split = log2(num_region); // number of dimension should be splited
	int combine_index = uniform(generator) * num_split; // randomly select an index in the splited dimension set
	int combine_dim = split_dim[combine_index]; // the dimension should be combined
	v1d<int> uncombine_region(num_region); // store the region which have not combine
	v1d<int> swallowed_region; // store the region which is swallowed up.

	// every region have not combine at first
	for (int j = 0; j < num_region; j++)
		uncombine_region[j] = j;

	// reduce the region until all region have combined
	while (uncombine_region.size() != 0)
	{
		// region 1 swallow region 2
		int region1 = uncombine_region.front();
		int region2 = region1 + pow(2, combine_index); // only one different dimension with region 1

		// reset the boundary of combine_dim dimension
		regions[region1].max_dim[combine_dim] = max_dim;
		regions[region1].min_dim[combine_dim] = min_dim;

		// combine the number of searchers
		regions[region1].num_searcher += regions[region2].num_searcher;

		// combine the sample solution
		samples[region1].insert(samples[region1].begin() + num_sample, 
								   samples[region2].begin(), samples[region2].begin() + num_sample);

		swallowed_region.push_back(region2); // store the region 2 which is swallowed by region 1
		uncombine_region.erase(uncombine_region.begin()); // region 1 have been combined
		uncombine_region.erase(std::find(uncombine_region.begin(), uncombine_region.end(), region2)); // region 2 have been combined
	}
	num_region = num_region / 2;
	num_sample = num_sample * 2;
	split_dim.erase(split_dim.begin() + combine_index);

	// 2. clear the reduced region
	for (int j = swallowed_region.size() - 1; j >= 0; j--)
	{
		int jj = swallowed_region[j];

		regions[jj].max_dim.clear();
		regions[jj].min_dim.clear();
		samples[jj].clear();

		regions.erase(regions.begin() + jj);
		samples.erase(samples.begin() + jj);

		expected_value.erase(expected_value.begin() + jj);
	}
}
void se::destroy()
{
	for (int j = 0; j < num_region; j++)
	{
		regions[j].max_dim.clear();
		regions[j].min_dim.clear();
		samples[j].clear();
	}
	regions.clear();
	samples.clear();

	expected_value.clear();

	searchers.clear();
	
	for (auto tmp : memory_param)
		tmp.clear();
	memory_param.clear();
	split_dim.clear();
}
