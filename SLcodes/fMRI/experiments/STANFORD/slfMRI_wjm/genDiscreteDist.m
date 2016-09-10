function [ discrete_dist ] = genDiscreteDist(pdf,n_trials)
	discrete_dist = round(pdf*n_trials)

end

