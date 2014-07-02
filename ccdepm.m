
function ccdepm(fname,cp,fnum)

	G = 0;
	S = [10,20,50];
	fe = 1000;	
	pop = cp.xmin + rand(cp.NP,cp.D) * (cp.xmax - cp.xmin);
	bestval = feval(fname,pop(1,:),fnum);
	[bestval,xbest] = evaluate(fname,cp,fnum,pop,bestval);
	while(G < cp.gen_max)
		

%%% Subgroup size choosing	
		tr1 = 0.6 - 0.5 * (G / cp.gen_max);
		tr2 = 0.9 - 0.5 * (G / cp.gen_max);
		for i = 1:cp.NP
			tr = rand(1);
			if(tr <= tr1)
				s = S(1);
			elseif(tr >= tr2)
				s = S(3);
			else 
				s = S(2);
			end
		end
%%%End Subgroup size choosing

		index = randperm(cp.D);
		for j = 1:(cp.D/s)
			l = (j-1)*s + 1;
			u = j * s;
			%fprintf(1,'Dimensions %d to %d\n',l,u);
			subpop = pop(:,index(l:u));
			subpop = adecm_for_ccdepm(fname,cp,fnum,xbest,bestval,subpop,index(l:u),fe,s);
			pop(:,index(l:u)) = subpop;
			[bestval,xbest] = evaluate(fname,cp,fnum,pop,bestval);
			
		end
%%% Perfunctory Mechanism
		r = zeros(1,cp.D);
		for i = 1:cp.NP
			pr = 0.05 + 0.05 * rand(1);
			if(rand <= pr)
				for q = 1:cp.D
					r(q) = floor(4 * rand(1));
				end
				bw = (cp.xmax - cp.xmin) * 0.1;
				for q = 1:cp.D
					if(r(q) == 0)
						pop(i,q) = pop(i,q) + bw * (-1 + rand * 2);
					end
				end
			end
		end
				
%%%End Perfunctory		
		[bestval,xbest] = evaluate(fname,cp,fnum,pop,bestval);
		fprintf(1,'GENERATION = %d, FITNESS VALUE = %d\n',G,bestval);
		G = G+1;					
	end
end

function [bestval,xbest] = evaluate(fname,cp,fnum,pop,bestval)
	xbest = pop(1,:);
	for i = 1:cp.NP
		xb = feval(fname,pop(i,:),fnum);
		if(xb <= bestval)
			xbest = pop(i,:);
			bestval = xb;
		end
	end
end		
