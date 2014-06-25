
function hs(fname,cp,fnum)

%%% INITIALIZATION
	HM = zeros(cp.HMS,cp.D);
	for i = 1:cp.HMS
		for j = 1:cp.D
			HM(i,j) = cp.xmin + rand(1)*(cp.xmax - cp.xmin);
		end
	end

%%%END INITIALIZATION

	impv = 0;
	HM_trial = zeros(cp.HMS,cp.D);
	selected = zeros(cp.HMS,cp.D);
	xb = feval(fname,HM(1,:),fnum);
	while( impv < cp.NI)
		
%%% IMPROVISATION
		for i = 1:cp.HMS
			for j = 1:cp.D
				if(rand(1) < cp.HMCR)
					temp_i = randi(cp.HMS);
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
		for i = 1:cp.HMS
			pv = feval(fname,HM(i,:),fnum);
			tv = feval(fname,HM_trial(i,:),fnum);
			if(pv  > tv)
				selected(i,:) = HM_trial(i,:);
            			se = tv;
			else 
				selected(i,:) = HM(i,:);
            			se = pv;
			end
		
			if(xb > se)
				xbest(1,:) = selected(i,:);
            			xb = se;
            end
                
        end
	
%%% END SELECTION
	impv = impv +1;
	HM = selected;


    %disp('generation='),disp(G),disp('smallest fitness value is'),disp(xb)
	fprintf(1,'fitness value = %d,no of improvisations = %d\n',xb,impv);
    end
end
	
