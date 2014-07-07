
function subpop = hs_for_subpop(fname,cp,fnum,xbest,bestval,subpop,index,fe,s,l,u)


	gpop = ones(cp.NP,1) * xbest;
	gpop(:,index) = subpop;
	F = zeros(1,cp.NP);
	fitness_vals = zeros(1,cp.NP);
	impv = 0;
	HM_trial = zeros(cp.NP,s);
	selected = zeros(cp.NP,s);
	G = 0;
	xb = bestval;
	
	while( G < 14)
		HM = subpop;
		G = G + 1;	
%%% IMPROVISATION
		for i = 1:cp.NP
			for j = 1:s
				if(rand(1) < cp.HMCR)
					temp_i = randi(cp.NP);
					HM_trial(i,j) = HM(temp_i,j);
					if(rand(1) < cp.PAR)
						if(rand(i) < 0.5)
							HM_trial(i,j) = HM_trial(i,j) + rand(1) * cp.bw;
						else 					
							HM_trial(i,j) = HM_trial(1,j) - rand(1) * cp.bw;
						end
					end
				else
				
					HM_trial(i,j) = cp.xmin + rand(1)*(cp.xmax - cp.xmin);
				end
			end
		end
		
		

%%%END IMPROVISATION		 


%%% SELECTION
		gpop_trial = gpop;
		gpop_trial(:,index) = HM_trial;
		pv = feval(fname,gpop,fnum);
		tv = feval(fname,gpop_trial,fnum);
		for i = 1:cp.NP
			if(pv(i)  > tv(i))
				selected(i,:) = HM_trial(i,:);
            			se = tv;
			else 
				selected(i,:) = HM(i,:);
            			se = pv;
			end
		
			if(xb > se)
            			xb = se;
            		end
                
        	end
	
%%% END SELECTION
	HM = selected;
	


    %disp('generation='),disp(G),disp('smallest fitness value is'),disp(xb)
	fprintf(1,'\n g = %g,fv = %d,dims = %g,from = %g,to = %g\n',G,xb,s,l,u);
	subpop = selected;
	end
end
	
