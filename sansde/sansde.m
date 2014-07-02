function sansde(fname,cp,fnum)
%     runs = 0;
%     xb_sum = 0;
% 	while(runs < 25)
%     
    G = 0;
    fe_count = 0;

%Generating Cauchy distribution for scale = 1 , location = 0    
    cauchy_arr = zeros(1,10000);
    for i = 1:10000
        cauchy_arr(i) = 1 / pi * ( 1 + i * i) ;
    end
    
%End Cauchy

    selected = zeros(cp.NP,cp.D);
    param_vectors = initialize_param_vectors(cp);
    xbest = param_vectors(1,:);
    xb = feval(fname,xbest,fnum);
    xb
    cp.fbias
	while((G < 1500) && (xb > cp.fbias))
		G = G + 1;
%%% MUTATION
        target_vectors = zeros(cp.NP,cp.D);
        for i = 1:cp.NP
            curr_val = feval(fname,param_vectors(i,:),fnum);
            fe_count = fe_count + 1;
            if(curr_val < xb)
                xb = curr_val;
                xbest = param_vectors(i,:);
            end
        end
        
        if(rem(G,50) == 1)
            ns1 = 0;  %mutation by rule1 and successful
            ns2 = 0;  %mutation by rule2 and successful
            nf1 = 0;  %mutation by rule1 and failed
            nf2 = 0;  %mutation by rule1 and failed
            
            frs1 = 0; % F chosen by rule 1 and successful
            frs2 = 0; % F chosen by rule 2 and successful
            frf1 = 0; % F chosen by rule 1 and failed
            frf2 = 0; % % F chosen by rule 2 and failed
            
            fp = 0.5;
            p = 0.5;
        else
            p = (ns1 * (ns2 + nf2)) / ((ns2 * (ns1 + nf1)) + (ns1 * (ns2 + nf2)));
            fp = (frs1 * (frs2 + frf2)) / ((frs2 * (frs1 + frf1)) + (frs1 * (frs2 + frf2)));
        end
        
        ns_rec = zeros(1,cp.NP);
        f_rec = zeros(1,cp.NP);
            
        
        for i = 1:cp.NP
            
%Calculating F
            if(rand(1) < fp)
                F = normrnd(0.5,0.3); % gaussion random number
                f_rec(i) = 1;
            else
                index = randi(10000,1);
                F = cauchy_arr(index);
                f_rec(i) = 2;
            end
% END calculating F
               
            r = generate_r(3,cp);
            xi = param_vectors(i,:);
            x_i1 = param_vectors(r(1),:);
            x_i2 = param_vectors(r(2),:);
            x_i3 = param_vectors(r(3),:);
            if(rand(1) < p )
                target_vectors(i,:) = x_i1 + F * (x_i2 - x_i3);
                ns_rec(i) = 1;
            else
                target_vectors(i,:) = xi + F *(xbest - xi) + F * (x_i1 - x_i2);
                ns_rec(i) = 2;
            end
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

%%% END MUTATION

%%% Crossover
		if(rem(G,25) == 1)
            if(G == 1)
                crm = 0.5;
            else
                crm = 0;
                for i = 1:count
                    crm = crm + w(i) * cr_rec(i);
                end
            end
            count = 0;
            cr_rec = zeros(1,25*cp.NP);
            w = zeros(1,25*cp.NP);
            diff_f = zeros(1,25*cp.NP);
            diff_f_sum = 0;
        end
        
         if(rem(G,5) == 1)
             cr = zeros(1,cp.NP);
         end
%             for i = 1:cp.NP
%                 cr_rec(cr_rec_count + i) = normrnd(crm,0.1);
%                 cr_sum = cr_sum + cr_rec(cr_rec_count + i);
%             end
%             cr_rec_count = cr_rec_count + cp.NP;
%         end
        
        trial_vectors = zeros(cp.NP,cp.D);
        
        for i = 1:cp.NP
            j_rand = round(1+ rand(1)*(cp.D - 1));
            if(rem(G,5) == 1)
                cr(i) = normrnd(crm,0.1);
            end
            for j = 1:cp.D
                if((rand(1) < cr(i)) ||  ( j == j_rand))
                    trial_vectors(i,j) = target_vectors(i,j);
                else
                    trial_vectors(i,j) = param_vectors(i,j);
                end
            end
        end
        
%%% END Crossover       
        for i = 1:cp.NP
            pv = feval(fname,param_vectors(i,:),fnum);
            tv = feval(fname,trial_vectors(i,:),fnum);
            if(pv  > tv)
                selected(i,:) = trial_vectors(i,:);
                if(rem(G,5) == 1)
                    count = count + 1;
                    diff_f(count) = pv - tv;
                    diff_f_sum = diff_f_sum + diff_f(count);
                    cr_rec(count) = cr(i);
                    w(count) = diff_f(count) / diff_f_sum;
                end
                
                if(ns_rec(i) == 1)
                    ns1 = ns1 + 1;
                else
                    ns2 = ns2 + 1;
                end
                if(f_rec(i) == 1)
                    frs1 = frs1 + 1;
                else
                    frs2 = frs2 + 1;
                end
                se = tv;
            else 
                selected(i,:) = param_vectors(i,:);
                if(ns_rec(i) == 1)
                    nf1 = nf1 + 1;
                else
                    nf2 = nf2 + 1;
                end
                if(f_rec(i) == 1)
                    frf1 = frf1 + 1;
                else
                    frf2 = frf2 + 1;
                end
                se = pv;
            end
		
            if(xb > se)
                xbest(1,:) = selected(i,:);
                xb = se;
            end
            fe_count = fe_count + 2;
                
        end
        param_vectors = selected;
        if(fe_count < cp.FE)
            fprintf(1,'generation = %d, fitness value = %d,fe = %d\n',G,xb,fe_count);
        end
    end
%     xb_sum = xb_sum + xb;
%     end
%     xb_mean = xb_sum /25
end
		
function param_vectors = initialize_param_vectors(cp)
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

