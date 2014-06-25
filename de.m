
function de(fname,cp,fn_num)
 
    format short
    G = 0;
    fe_count = 0;
    param_vectors = generate_vectors( cp);
    selected = zeros(cp.NP,cp.D);
    xbest = param_vectors(1,:);
    xb = feval(fname,xbest,fn_num);
    while(1)    
        G = G + 1;
        target_vectors = DE_rand_1(param_vectors,cp);
        trial_vectors = exponential_crossover(param_vectors,target_vectors,cp);
        %trial_vectors = binomial_crossover(param_vectors,target_vectors,cp);
% [xb,xbest,param_vectors] = selection(fname,param_vectors,trial_vectors,cp);
	for i = 1:cp.NP
		pv = feval(fname,param_vectors(i,:),fn_num);
		tv = feval(fname,trial_vectors(i,:),fn_num);
		if(pv  > tv)
			selected(i,:) = trial_vectors(i,:);
            se = tv;
		else 
			selected(i,:) = param_vectors(i,:);
            se = pv;
		end
		
		if(xb > se)
			xbest(1,:) = selected(i,:);
            xb = se;
		end
        fe_count = fe_count + 2;
                
    end
        param_vectors = selected;
        fprintf(1,'generation = %d, fitness value = %g,fe = %d\n',G,xb,fe_count);
	
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

function r = generate_r(num,cp)
	mu = 0;
	r = zeros(1,num);
	while(mu == 0)
        for i = 1:num
			r(i) = round(1 + rand(1) * (cp.NP - 1));
        end
        if (size(r,2) == size(unique(r),2))
			mu = 1;
        end
	end

end

function target_vectors = DE_rand_1(param_vectors,cp)
	target_vectors = zeros(cp.NP,cp.D);
	for i = 1:cp.NP
		r = generate_r(3,cp);
		r1 = param_vectors(r(1),:);
		r2 = param_vectors(r(2),:);
		r3 = param_vectors(r(3),:);
		target_vectors(i,:) = r1 + cp.F * (r2 - r3);
	end
	for i = 1: cp.NP
		for j = 1: cp.D
            if(target_vectors(i,j) < cp.xmin)
				target_vectors(i,j) = cp.xmin;
            end
			if(target_vectors(i,j) > cp.xmax)
				target_vectors(i,j) = cp.xmax;
			end	
		end
	end

end    

function L = calculate_l(cp)

	L = 1;
    while((rand(1) < cp.Cr) && (L < cp.D))
		L = L + 1;
    end
end


function trial_vectors = exponential_crossover(param_vectors,target_vectors,cp)

	trial_vectors = zeros(cp.NP,cp.D);
	for i = 1:cp.NP
		n = round(1+ rand(1)*(cp.D - 1));
		L = calculate_l(cp);
		for j = 1:cp.D
			if ((j >= n) && (j < n + L))
				trial_vectors(i,j) = target_vectors(i,j);
			else
				trial_vectors(i,j) = param_vectors(i,j);
			end
		end
	end
end

function trial_vectors = binomial_crossover(param_vectors,target_vectors,cp)

	trial_vectors = zeros(cp.NP,cp.D);
	for i = 1:cp.NP
		j_rand = round(1+ rand(1)*(cp.D - 1));
		for j = 1:cp.D
			if((rand(1) < cp.Cr) ||  ( j == j_rand))
				trial_vectors(i,j) = target_vectors(i,j);
			else
				trial_vectors(i,j) = param_vectors(i,j);
			end
		end
	end
end


% function [xb,xbest,selected] = selection(fname,param_vectors,trial_vectors,cp)
% 	
% 	selected = zeros(cp.NP,cp.D);
% 	xbest = param_vectors(1,:);
% 	for i = 1:cp.NP
% 		pv = feval(fname,param_vectors(i,:));
% 		tv = feval(fname,trial_vectors(i,:));
% 		if(pv  > tv)
% 			selected(i,:) = trial_vectors(i,:);
% 		else
% 			selected(i,:) = param_vectors(i,:);
% 		end
% 		xb = feval(fname,xbest);
% 		se = feval(fname,selected(i,:));
% 		if(xb > se)
% 			xbest(1,:) = selected(i,:);
% 		end
% 	end
% end
