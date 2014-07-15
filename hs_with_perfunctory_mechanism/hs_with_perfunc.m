
function bestval = hs_with_perfunc(fname,cp,fnum,run)

	G = 0;
	S = [10,20,50];
	fe = 1000;	
	pop = cp.xmin + rand(cp.NP,cp.D) * (cp.xmax - cp.xmin);
	bestval = feval(fname,pop(1,:),fnum);
	[bestval,xbest] = evaluate(fname,cp,fnum,pop,bestval);
	filename = strcat(int2str(fnum),'_',int2str(run));
	f = fopen(filename,'w');
	while(G < cp.gen_max)
		
%%% Subgroup size choosing	
%		tr1 = 0.6 - 0.5 * (G / cp.gen_max);
%		tr2 = 0.9 - 0.5 * (G / cp.gen_max);
%		for i = 1:cp.NP
%			tr = rand(1);
%			if(tr <= tr1)
%				s = S(1);
%			elseif(tr >= tr2)
%				s = S(3);
%			else 
%				s = S(2);
%			end
%		end
		s=10;
%%%End Subgroup size choosing

		index = randperm(cp.D);
		for j = 1:(cp.D/s)
			l = (j-1)*s + 1;
			u = j * s;
			%fprintf(1,'Dimensions %d to %d\n',l,u);
			k = index(l:u);
			subpop = pop(:,k);
			subpop = hs_for_subpop(fname,cp,fnum,xbest,bestval,subpop,k,fe,s,l,u);
			pop(:,index(l:u)) = subpop;
			[bestval,xbest,ind] = evaluate(fname,cp,fnum,pop,bestval);
			
		end
		bestval
%%% Perfunctory Mechanism
		r = zeros(cp.NP,cp.D);
%		for i = 1:cp.NP
%			pr = 0.05 + 0.05 * rand(1);
%			if(rand(1) <= pr)
%				for q = 1:cp.D
%					r(q) = floor(4 * rand(1));
%				end
%				bw = (cp.xmax - cp.xmin) * 0.1;
%				for q = 1:cp.D
%					if(r(q) == 0)
%						pop(i,q) = pop(i,q) + bw * (-1 + rand * 2);
%					end
%				end
%			end
%		end
pr = 0.05 + 0.05 * rand(cp.NP,1);
pr_rep = repmat(pr,1,cp.D);
rand_temp = rand(cp.NP,1);
rand_temp_rep = repmat(rand_temp,1,cp.D);
p_true = rand_temp_rep <= pr_rep;
p_false = rand_temp_rep > pr_rep;
%size(p_true)
%size(floor(4*rand(cp.NP,cp.D)))
r = p_true .*  floor(4 * rand(cp.NP,cp.D)) ;
r_zero = (r==0);
bw = ones(cp.NP,cp.D) * (cp.xmax - cp.xmin) * 0.1;
pop = p_true .* r_zero .* (pop + bw * (-1 + rand * 2)) + p_false .* pop;

perfunc_best = feval(fname,pop(ind,:),fnum);
if(perfunc_best > bestval)
	pop(ind,:) = xbest;
end
 				
%%%End Perfunctory		
		[bestval,xbest,index] = evaluate(fname,cp,fnum,pop,bestval);
		fprintf(f,'GENERATION = %d, FITNESS VALUE = %d\n',G,bestval);
%		fprintf(1,'\n GENERATION = %d, FITNESS VALUE = %d\n',G,bestval);
		prev_bestval = bestval;

		G = G+1;					
	end
	fclose(f);
end

function [bestval,xbest,index] = evaluate(fname,cp,fnum,pop,bestval)
		xb = feval(fname,pop,fnum);
		[bestval,index] = min(xb);
		xbest = pop(index,:);
end		
