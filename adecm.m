
function adecm(fname,cp,fnum)

	fl = 0.5;
	fs = 0.15;
	cr = 0.9;
	param_vectors = generate_vectors(cp)
	F = zeros(1,cp.NP);
	fitness_vals = zeros(1,cp.NP);
	for i = 1:cp.NP
		fitness_vals(i) = feval(fname,param_vectors(i,:),fnum);
	end
	for G = 1:cp.gen_max
		sf = [];
		C = 0.8 - 0.6 * (G/cp.gen_max);
		[fitness_vals,param_vectors] = sort_by_fitness(fitness_vals,param_vectors);
		for i = 1:cp.NP
			F(i) = rand_cauchy(fl,fs);
			[r1,r2,r3] = generate_r_selective_pressure(param_vectors,B,cp,i);
			x_r1 = param_vectors(r1,:);
			x_r2 = param_vectors(r2,:);
			x_r3 = param_vectors(r3,:);
			xbest = param_vectors(1,:);
			if(C < rand(1))
				target_vectors(i,:) = x_r1 + F(i) * (x_r2 - x_r3);
			else
				target_vectors(i,:) = x_r1 + F(i) * (xbest - param_vectors(i,:) +  F(i) * (x_r2 - x_r3);
			end
			j_rand = randi(cp.D);
			for j = 1:cp.D
				if((j==j_rand) || (rand(1) <= cr) 
					trial_vectors(i,j) = target_vectors(i,j);
				else
					trial_vectors(i,j) = param_vectors(i,j);
				end
			end
			tv = feval(fname,trial_vectors(i,:))
			if(tv <= fitness_vals(i))
				param_vectors(i,:) = trial_vectors(i,:);
				fitness_vals(i) = tv;
				sf(end+1) = F(i);
			end
		end
				
	end
end


function param_vectors = generate_vectors(cp)
	param_vectors = zeros(cp.NP,cp.D);
	for i = 1: cp.NP
		for j = 1: cp.D
			param_vectors(i,j) = cp.xmin + rand(1)*(cp.xmax - cp.xmin);
		end
	end
end

function [fitness_vals,param_vectors] = sort_by_fitness(fitness_vals,param_vectors)
	
	fitness_params = [fitness_vals param_vectors];
	fitness_params = sortrows(fitness_params);
	param_vectors = fitness_params(:,2:end);
	fitness_vals = fitness_params(:,1);
end	


function [r1,r2,r3] = generate_r_selective_pressure(param_vectors,B,cp,i)

	r1 = randi(cp.NP);
	r2 = randi(cp.NP);
	while((r1 == r2)||(r2 = i) || (r1 == i))
		r2 = randi(cp.NP);
                r1 = randi(cp.NP);
	end
	r3 = randi([floor(cp.NP * 0.8),cp.NP]);
	while((r3 == r1)|| (r3 == r2) || (r3 == i))
		r3 = randi([floor(cp.NP * 0.8),cp.NP]);
	end	
end
