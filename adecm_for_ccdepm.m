
function param_vectors = adecm_for_ccdepm(fname,cp,fnum,xbest,bestval,param_vectors,index,fe,s)

	gpop = ones(cp.NP,1) * xbest;
	gpop(:,index) = param_vectors;
	target_vectors = zeros(cp.NP,s); 
	trial_vectors = zeros(cp.NP,s); 
	F = zeros(1,cp.NP);
	fitness_vals = zeros(1,cp.NP);
	for i = 1:cp.NP
		fitness_vals(i) = feval(fname,gpop(i,:),fnum);
	end
	for G = 1:30
		s_f = zeros(1,cp.NP);
        	s_f_count = 0;
		C = 0.8 - 0.6 * (G/cp.gen_max);
		[fitness_vals,param_vectors] = sort_by_fitness(fitness_vals,param_vectors);
		bestval_now = fitness_vals(1);
		for i = 1:cp.NP
			F(i) = abs(rand_cauchy(cp.fl,cp.fs));
			[r1,r2,r3] = generate_r_selective_pressure(cp,i);
			x_r1 = param_vectors(r1,:);
			x_r2 = param_vectors(r2,:);
			x_r3 = param_vectors(r3,:);
			xbest = param_vectors(1,:);
			if(C < rand(1))
				target_vectors(i,:) = x_r1 + F(i) * (x_r2 - x_r3);
			else
				target_vectors(i,:) = x_r1 + F(i) * (xbest - param_vectors(i,:)) +  F(i) * (x_r2 - x_r3);
			end
			j_rand = randi(s);
			for j = 1:s
				if((j==j_rand) || (rand(1) <= cp.cr)) 
					trial_vectors(i,j) = target_vectors(i,j);
				else
					trial_vectors(i,j) = param_vectors(i,j);
				end
			end
			g_trial = ones(cp.NP,1) * xbest;
			g_trial(:,index) = trial_vectors;
			tv = feval(fname,g_trial(i,:),fnum);
			if(tv <= fitness_vals(i))
				param_vectors(i,:) = trial_vectors(i,:);
				fitness_vals(i) = tv;
                		s_f_count = s_f_count + 1;
				s_f(s_f_count) = F(i);
				if(tv < bestval)
					bestval = tv;
				end
			end
		end
		gpop(:,index) = param_vectors;
		pmean_sf = nthroot(sum(s_f.^1.5)/length(s_f),1.5);
		cf = 0.9 + rand(1) * 0.1;
		cp.fl = cf * cp.fl + (1-cf)* pmean_sf;
		cp.fs = cf * cp.fs + (1-cf)* pmean_sf;
%        	fprintf(1,'generation = %d, fitness value = %d,error = %d\n',G,bestval,cp.fbias - bestval);
	end
				
end



function [fitness_vals,param_vectors] = sort_by_fitness(fitness_vals,param_vectors)
	
	fitness_params = [fitness_vals' param_vectors];
	fitness_params = sortrows(fitness_params);
	param_vectors = fitness_params(:,2:end);
	fitness_vals = fitness_params(:,1)';
end	


function [r1,r2,r3] = generate_r_selective_pressure(cp,i)

	r1 = floor(rem(cp.NP/(2*(cp.B - 1))* sqrt(cp.B * cp.B - 4 * (cp.B -1)*rand(1)),cp.NP));
	r2 = floor(rem(cp.NP/(2*(cp.B - 1))* sqrt(cp.B * cp.B - 4 * (cp.B -1)*rand(1)),cp.NP));

	while((r1 == r2)||(r2 == i) || (r1 == i)|| (r1*r2 == 0))
		r1 = floor(rem(cp.NP/(2*(cp.B - 1))* sqrt(cp.B * cp.B - 4 * (cp.B -1)*rand(1)),cp.NP));
		r2 = floor(rem(cp.NP/(2*(cp.B - 1))* sqrt(cp.B * cp.B - 4 * (cp.B -1)*rand(1)),cp.NP));
	end
	r3 = randi([floor(cp.NP * 0.8),cp.NP]);
	while((r3 == r1)|| (r3 == r2) || (r3 == i)|| (r3 == 0))
		r3 = randi([floor(cp.NP * 0.8),cp.NP]);
	end	
end

function num = rand_cauchy(location, scale)

	p = 0.0;
	while (p == 0.0)
		p = rand(1);
		num = location + scale *  tan( pi * (p - 0.5));
	end
end
