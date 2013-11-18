%%%%%%%%%%%%%%%%%%%%% TEK292: LAB - Genetic Algorithms %%%%%%%%%%%%%%%%%%%%%%%%
%%
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
%% options		- [selection_flag_p2, recombination_flag, proximity_flag, plot_flag]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function: GA.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [top_ind] = GA(pop, n_gen,  fit_func,  func_interval, ...
		crossover_prob, mutate_rate, elitism_rate, options)

	%%%%%%%%% Algorithm stucture %%%%%%%%%
	% - initiate
	% - set options 
	% - Set selection function
	% - Set recombination function

	% - Start evolution
	%	- Get top quantile
	%		- Selection
	%		- Recombination
	%	- Mutate
	%	- Recalculate fitness
	%	- Plot	
	%	- Check convergence
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if nargin < 7
		error('input_example :  first 7 inputs are required')
	end
	if nargin < 8
		options = []
	end

	% Initiate
	% set options
	[selection_flag_p2, recombination_flag, proximity_flag, plot_flag] = set_options(options);
	[pop_size n]	=	 size(pop);
	pop 		=	 sortrows(pop, [2]);
	n_genes		=	 n-3;
	n_spicies	=	 max(pop(:,end));
	top_quantile	=	 ceil(pop_size*0.25);
	top_fitness	=	 cell(1, n_spicies);
	avg_fitness	=	 cell(1, n_spicies);
	top_individ	=	 cell(1, n_spicies);
	conv_flag	=	 0
	
	for i=1:1:n_spicies
		s_pop = spiecies_pop(pop, [0 0 0 0 i]);
		tot_fitness{1, i} = [sum(s_pop(:,2))];
		avg_fitness{1, i} = [mean(s_pop(:,2))];
		top_individ{1, i} = [pop(1, 2)];
	end
	tot_fitness{1,1}
	% Set up plotting	
	if plot_flag >= 1  	
		[X, Y] 	= meshgrid(func_interval(1):0.2:func_interval(2));
		Z 	= fit_func(X,Y);
	end	         	

	% Set selection function for parent 1
	switch selection_flag_p2
		case 1 % random
			selection_func_p2 = @rand_ind;
		case 2 % Probability wheel
			selection_func_p2 = @prob_wheel;
		case 3 % Rank based
			selection_func_p2 = @rank_based;

	end

	% Set recombination function
	switch recombination_flag
		case 1 % Average
			recombination_func = @avg_cross;
		case 2 % Extrapolation
			recombination_func = @extrapol_cross;
	end
	
	% TODO
	% Set insertion function
	% switch insertion_flag
	%	case 1 % From bottom
	%
	%	case 2 % keep bottom 10%

	% selection_func_p1
	selection_func_p2
	recombination_func
	proximity_flag
	plot_flag
	
	%% Start evolution
	for gen=1:1:n_gen
		% Selection & recombination for top individuals
		if rand(1) < crossover_prob
			for i=1:1:top_quantile
				%% Selection
				% p1 - parent 1
				% p2 - parent 2
				p1 = pop(i, :);
				% Remove p1 from selection population to prevent
				% inbreeding
				
				new_pop = spiecies_pop(remove_ind(pop, p1), p1);
				% If selection is proximity based
				if proximity_flag==1
					new_pop = near_pop(pop, p1);
				end
				p2 = selection_func_p2(new_pop);

				%% recombination
				% o - offspring genes
				o = recombination_func(p1, p2, fit_func, func_interval);
				% insert new individual at the bottom of the population
				pop(pop_size-i+1, 3:4) = o(:);
			end
		end
		
		% Mutation
		% We do not wish to mutate the top individuals (elitism)
		% n_mutate - genes to mutate 
		n_mutate = ceil(mutate_rate*pop_size*n_genes);
		for i=1:1:n_mutate
			ind = randi([elitism_rate pop_size], 1);
			gene = randi([3 (2+n_genes)], 1);
			if rand(1)<0.5
				pop(ind, gene) = pop(ind, gene)+0.1; 
			else
				pop(ind, gene) = pop(ind, gene)-0.1;
			end
		end	

		% Recalculate fitness & sort population
		pop(:,2) = fit_func(pop(:,3), pop(:,4));
		pop = sortrows(pop, [2]);
		
	
		for s=1:1:n_spicies	
			s_pop = spiecies_pop(pop, [0 0 0 0 s]);
			tot_fitness{1, s} = [tot_fitness{1, s} sum(s_pop(:, 2))];
			avg_fitness{1, s} = [avg_fitness{1, s} mean(s_pop(:, 2))];
			top_individ{1, s} = [top_individ{1, s} s_pop(1, 2)];
		end
		top_individ = [top_individ pop(1, 2)];
		top_ind	    = pop(1,:);
		% plot
		if plot_flag >= 1
			figure(1);
			subplot(211)
			mesh(X,Y,Z);
			hold on
			col = 'bgrcmyk';
			for i=1:1:pop_size
				%plotc = 'r*';
				%[pop(i,3), pop(i,4), pop(i,2)]
				%[col(pop(i,end))]
				plot3(pop(i,3), pop(i,4), pop(i,2),  [col(pop(i,end)), '*']);
			end
			hold off
			ylabel(['Generation: ', num2str(gen)]);
			subplot(212)
			hold on
			for s=1:1:n_spicies
				%top_individ{1, s}
				plot([1:1:gen+1], avg_fitness{1, s}, [col(s), '-']);
				plot([1:1:gen+1], top_individ{1, s}, [col(s), '--']);
			end
			legend('Average fitness', 'Top individ');
			hold off
			
			%if conv_flag==0
			%	xlabel(['Press key to continue']);
			%	pause;
			%end

		end

		% Check convergence
		tol = 0.05;
		if abs(pop(1, 2) - mean(pop(:, 2))) < tol
			conv_gen = gen;
			conv_flag = 1;	
		end
		if conv_flag == 1
			break;
		end
	end
end


% TODO
% - Look over proximity selection
% - Quick convergence for top individual, then hard to find other minimum
% - 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Help functions      %%%%%%%%%%%%%%%%%%%%%%%%%%%%
function new_pop = remove_ind(pop, ind)
	new_pop = pop((pop(:,1)==ind(1))==0, :); 
end
% Get population with same spieces
function new_pop = spiecies_pop(pop, ind)
	new_pop = pop((pop(:,end)==ind(end)), :);	
end
% 
function new_pop = near_pop(pop, ind)

	[len dim] = size(pop);

	D = NaN(len-1, 2);
	for i=2:1:len
		D(i-1, :) = [pop(i,1), sqrt(sum(((ind(3:4)-pop(i,3:4)).^2)))];
	end
	D_sorted = sortrows(D, [2]);
	thresh = ceil(len*0.25);
	new_pop = pop(D_sorted(1:thresh, 1), :);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Selection functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% arguments:
% pop 		- a population 
%% returns:
% new_pop 	- the population without the selected individual
% ind		- selected individual

function [ind, new_pop] = rank_based(pop)
    ind = pop(1, :)	
    new_pop = pop((pop(:,1)==ind(1))==0, :); 
end

function [ind, new_pop] = rand_ind(pop)
	[m n] = size(pop);
	ind = pop(randi([1 m], 1), :);
	new_pop = pop((pop(:,1)==ind(1))==0, :); 
end

function [ind, new_pop] = prob_wheel(pop)

	pop = sortrows(pop, [2]);
	s = sum(pop(:,2));

	prob = [pop(:,2)./s];

	c = cumsum(prob);

	i = 0;
	r = rand(1);
	temp = find((c < r));
	if isempty(temp)
		i = 1;
	else
		i = temp(end)+1;
	end
	ind = pop(i, :);
	new_pop = pop((pop(:,1)==ind(1))==0, :); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Recombination functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% arguments:
% p1 		- parent 1
% p2 		- parent 2 
%% returns:
% offspring 	- genes for offspring of p1 & p2

function offspring = avg_cross(p1, p2, ~, ~)
	avg = @(x1, x2) (x1 + x2)/2;
	offspring = [avg(p1(3), p2(3)) avg(p1(4), p2(4))];
end

function offspring = extrapol_cross(p1, p2, fitfunc, interval)
	c = 0.5;
	o1 = [c*p1(3) + c*p2(3), 3*c*p1(3) - c*p2(3), -c*p1(3) + 3*c*p2(3)];
	o2 = [c*p1(4) + c*p2(4), 3*c*p1(4) - c*p2(4), -c*p1(4) + 3*c*p2(4)];
	index1 = find(abs(o1)<interval(2));
	index2 = find(abs(o2)<interval(2));
	o1 = o1(index1);
	o2 = o2(index2);

	fitness = [];

	for i=1:1:length(o1)
		for k=1:1:length(o2)
			fitness = [fitness; i k fitfunc(o1(i), o2(k))];
		end
	end
	% Minimize (2) / Maximize (2)
	fitness = sortrows(fitness, [2]);
	offspring = [o1(fitness(1, 1)) o2(fitness(1,2))];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [selection_flag_p2, recombination_flag, proximity_flag, plot_flag] = set_options(options)
	% Standard values:
	selection_flag_p2	= 1;
	recombination_flag 	= 1;
	proximity_flag		= 0;
	plot_flag	 	= 0; 

	for i=1:1:length(options)
		if i == 1
			selection_flag_p2 = options(1);
		elseif i == 2
			recombination_flag = options(2);
		elseif i == 3
			proximity_flag = options(3);
		elseif i == 4
			plot_flag = options(4);
		end
	end
end
