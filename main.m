%%%%%%%%%%%%%%%%%%%%% TEK292: LAB - Genetic Algorithms %%%%%%%%%%%%%%%%%%%%%%%%
% expertmaker.org/tek292
% lars@expertmaker.com
%% by: Ludwig Jacobsson | knd09lja | ludjac@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters and keywords:
%% pop 			- population
%% n_gen		- number of generation
%% fit_func		- fitness function
%% func_interval	- interval for function to optimize (used for plotting)
%% crossover_prob 	- Probability for mating
%% mutate_rate		- Rate for mutation
%% elitism_rate		- Rate for elitism
%% options		- [selection_flag, recombination_flag, spicies_flag,  plot_flag]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function: main.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% files: 
% main.m
% GA.m
%% 
%%%%%%%%% Algorithm stucture %%%%%%%%%
% initiate
	% set options 
% Set recombination function
% Set selection function

% Start evolution
	% Get top quantile
		% Selection
		% Recombination
	% Mutate
	% Recalculate fitness
	% Plot	
	% Check convergence

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initiate 
pop_size        = 200;
n_genes		= 2;
n_gen		= 100;
%fitfunc         = @(x, y) 1./(1 + x.^2 + y.^2);
%interval        = [-3 3];
fit_func        = @(x, y) y.*sin(4*x)+1.1*x.*sin(2*y);
func_interval   = [0 10];
mutate_rate     = 0.01;
elitism_rate    = 3;
crossover_prob  = 0.9;
n_species       = 3;

%%% Create population with matrix layout as seen bellow

% pop= | indN  | Fitness      | gene1         | gene2          (| geneN )                  | species                       |

pop = [[1:pop_size]; zeros(1,pop_size); ...
	[func_interval(2) + (func_interval(1)-func_interval(2)).*rand(n_genes,pop_size)] ...
	;randi([1 n_species], 1, pop_size) ]';
% pop = [[1:Npop]; zeros(1,Npop); [interval(2) + (interval(1)-interval(2)).*rand(Ngen,Npop)] ;randi([1 Nspecies], 1, Npop) ]';
% Generate fitness for population
pop(1:end, 2) = fit_func(pop(1:end, 4), pop(1:end, 5));

%function [avg,  top, conv_gen] = GA(pop, n_gen,  fit_func,  func_interval, ...
%	crossover_prob, mutate_rate, elitism_rate, options)
%% options: 
% [selection_flag_p2, recombination_flag, proximity_flag, plot_flag]
% Selection:
% 1 - random
% 2 - Probability wheel
% 3 - Rank based
% Recombination
% 1 - Average
% 2 - Extrapolation
% Proximity
% 1 - Selection is based on who is nearby
% 0 - proximity OFF
% Plot flag:
% 1 - plots
% 0 - no plots
options = [2, 1, 1, 1];
[top_ind] = GA(pop, n_gen, fit_func, func_interval, crossover_prob, mutate_rate, elitism_rate, options)
