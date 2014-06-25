
    % F_VTR		"Value To Reach" (stop when ofunc < F_VTR)
	F_VTR =120;

% I_D		number of parameters of the objective function 
		I_D = 10; 

% FVr_minbound,FVr_maxbound   vector of lower and bounds of initial population
%    		the algorithm seems to work especially well if [FVr_minbound,FVr_maxbound] 
%    		covers the region where the global minimum is expected
%               *** note: these are no bound constraints!! ***
      FVr_minbound = -5*ones(1,I_D); 
      FVr_maxbound = 5*ones(1,I_D);     
      d=1;
            
% I_NP            number of population members
		I_NP = 100;  %pretty high number - needed for demo purposes only

% I_itermax       maximum number of iterations (generations)
		I_itermax =1005; 
       
% F_weight        DE-stepsize F_weight ex [0, 2]
		 

% F_CR            crossover probabililty constant ex [0, 1]
		F_CRMAX=1;
        F_CRMIN=0.5;
        
        F_CR=0.5;
     code=input('plz enter the no. of the function u want to evaluate:->');
 
% I_strategy        
     I_strategy = input('Enter which DE mutation scheme to use(crossover is binomial type)\n\n1 --> DE/rand/1\n\n2 --> DE/best/1\n\n 3 --> DE/target-to-best/1\n\n4 --> DE/best/2\n\n 5 --> DE/rand/2\n\n 6 --> DE/target-to-best/p\n\n 7 --> DE/trigonmrtric-mutation\n\n:');
     if (I_strategy==6)
         P=input('ENTER P LYING BETWEEN 0 AND 100 -->');
     end
% I_refresh     intermediate output will be produced after "I_refresh"
%               iterations. No intermediate output will be produced
%               if I_refresh is < 1
      I_refresh = input('After how many iterations you want to see the output:');  
           
      %-----Check input variables--------------------
if (I_NP < 5)
   I_NP=5;
   fprintf(1,' I_NP increased to minimal value 5\n');
end
if ((F_CR < 0) | (F_CR > 1))
   F_CR=0.5;
   fprintf(1,'F_CR should be from interval [0,1]; set to default value 0.5\n');
end
if (I_itermax <= 0)
   I_itermax = 200;
   fprintf(1,'I_itermax should be > 0; set to default value 200\n');
end
I_refresh = floor(I_refresh);

%-----Initialize population and some arrays-------------------------------
FM_pop = zeros(I_NP,I_D); %initialize FM_pop to gain speed

%----FM_pop is a matrix of size I_NPx(I_D+1). It will be initialized------
%----with random values between the min and max values of the-------------
%----parameters-----------------------------------------------------------

for k=1:I_NP
   FM_pop(k,:) = FVr_minbound + rand(1,I_D).*(FVr_maxbound - FVr_minbound);
end

FM_popold     = zeros(size(FM_pop));  % toggle population
FVr_bestmem   = zeros(1,I_D);% best population member ever
FVr_bestmemit = zeros(1,I_D);% best population member in iteration
I_nfeval      = 0;                    % number of function evaluations

%------Evaluate the best member after initialization----------------------

I_best_index   = 1;                   % start with first population member
F_val(1)       = benchmark_func(FM_pop(I_best_index,:),code);

F_bestval = F_val(1);                 % best objective function value so far
I_nfeval  = I_nfeval + 1;

for k=2:I_NP                          % check the remaining members
  F_val(k)  = benchmark_func(FM_pop(k,:),code);
  I_nfeval  = I_nfeval + 1;
  if (left_win(F_val(k),F_bestval) == 1)
     I_best_index   = k;              % save its location
     F_bestval      = F_val(k);
  end   
end
FVr_bestmemit = FM_pop(I_best_index,:); % best member of current iteration
F_bestvalit   = F_bestval;              % best value of current iteration
F_maxerr=F_bestval-F_VTR;

FVr_bestmem = FVr_bestmemit;            % best member ever

%FM_normal=0.4+.3*randn(1200,1);
%FM_sortednormal1=sort(FM_normal,'descend');
%i=1;
%for h=1:1200
 %if(FM_sortednormal1(h,1)>=0 & FM_sortednormal1(h,1)<0.94)
   %    FM_sortednormal(i,1)=FM_sortednormal1(h,1);
  %     i=i+1;
 %   end
%end


%------DE-Minimization---------------------------------------------
%------FM_popold is the population which has to compete. It is--------
%------static through one iteration. FM_pop is the newly--------------
%------emerging population.----------------------------------------

FM_sortedpopold=zeros(I_NP,I_D);
FM_pm1   = zeros(I_NP,I_D);   % initialize population matrix 1
FM_pm2   = zeros(I_NP,I_D);   % initialize population matrix 2
FM_pm3   = zeros(I_NP,I_D);   % initialize population matrix 3
FM_pm4   = zeros(I_NP,I_D);   % initialize population matrix 4
FM_pm5   = zeros(I_NP,I_D);   % initialize population matrix 5
FM_bm    = zeros(I_NP,I_D);   % initialize FVr_bestmember  matrix
FM_ui    = zeros(I_NP,I_D);   % intermediate population of perturbed vectors
FM_mui   = zeros(I_NP,I_D);   % mask for intermediate population
FM_mpo   = zeros(I_NP,I_D);   % mask for old population
FVr_rot  = (0:1:I_NP-1);               % rotating index array (size I_NP)
FVr_rotd = (0:1:I_D-1);       % rotating index array (size I_D)
FVr_rt   = zeros(I_NP);                % another rotating index array
FVr_rtd  = zeros(I_D);                 % rotating index array for exponential crossover
FVr_a1   = zeros(I_NP);                % index array
FVr_a2   = zeros(I_NP);                % index array
FVr_a3   = zeros(I_NP);                % index array
FVr_a4   = zeros(I_NP);                % index array
FVr_a5   = zeros(I_NP);                % index array
FVr_ind  = zeros(4);

FM_meanv = ones(I_NP,I_D);

I_iter = 1;
p=F_bestval;
while ((I_iter<I_itermax) & (F_bestval> F_VTR+1.e-15 | F_bestval< F_VTR-1.e-15))
    
    %p=100-0.099*I_iter;
    p=25;
    
  FM_popold = FM_pop;  % save the old population
  F_valsort=sort(F_val);
  for h=1:I_NP
  for g=1:I_NP
      if(F_val(g) == F_valsort(h))
          FM_sortedpopold(h,:)=FM_popold(g,:);
          break;
      end
  end
  end
  FVr_ind = randperm(4);               % index pointer array

  FVr_a1  = randperm(I_NP);                   % shuffle locations of vectors
  FVr_rt  = rem(FVr_rot+FVr_ind(1),I_NP);     % rotate indices by ind(1) positions
  FVr_a2  = FVr_a1(FVr_rt+1);                 % rotate vector locations
  FVr_rt  = rem(FVr_rot+FVr_ind(2),I_NP);
  FVr_a3  = FVr_a2(FVr_rt+1);                
  FVr_rt  = rem(FVr_rot+FVr_ind(3),I_NP);
  FVr_a4  = FVr_a3(FVr_rt+1);               
  FVr_rt  = rem(FVr_rot+FVr_ind(4),I_NP);
  FVr_a5  = FVr_a4(FVr_rt+1);                

  FM_pm1 = FM_popold(FVr_a1,:);             % shuffled population 1
  FM_pm2 = FM_popold(FVr_a2,:);             % shuffled population 2
  FM_pm3 = FM_popold(FVr_a3,:);             % shuffled population 3
  FM_pm4 = FM_popold(FVr_a4,:);             % shuffled population 4
  FM_pm5 = FM_popold(FVr_a5,:);             % shuffled population 5

  for k=1:I_NP                              % population filled with the best member
    FM_bm(k,:) = FVr_bestmemit;             % of the last iteration
  end
  
  
  
 FM_mui = rand(I_NP,I_D) < F_CR;  % all random numbers < F_CR are 1, 0 otherwise
 FM_mpo = FM_mui < 0.5;    % inverse mask to FM_mui

  if (I_strategy == 1)    % DE/rand/1
           F_weight=0.8;
           FM_ui = FM_pm3 + F_weight*(FM_pm1 - FM_pm2);   % differential variation
           FM_ui = FM_popold.*FM_mpo + FM_ui.*FM_mui;     % crossover
           FM_origin = FM_pm3;
  elseif (I_strategy == 2)% DE/target-to-best/1
           for k=1:I_NP 
           for j=1:I_D    
           F_weight=0.5+0.3*randn;
                  FM_ui(k,j) = FM_popold(k,j) +0.85*(FM_bm(k,j)-FM_popold(k,j)) +F_weight*(FM_pm1(k,j) - FM_pm2(k,j));
           end      
           q=ceil(p*rand(1,1));
           for j=1:I_D
                   if(rand(1,1)>F_CR)
                       FM_ui(k,j)=FM_sortedpopold(q,j);
                   end
           end
           end
  elseif (I_strategy == 3)                         % DE/best/1 
            for k=1:I_NP 
           for j=1:I_D
                  F_weight=0.5+0.3*randn;
                  FM_ui(k,j) = FM_bm(k,j) +F_weight*(FM_pm1(k,j) - FM_pm2(k,j));
           end
           q=ceil(p*rand(1,1));
           for j=1:I_D
                   if(rand(1,1)>F_CR)
                       FM_ui(k,j)=FM_sortedpopold(q,j);
                   end
           end
           end
         
  elseif (I_strategy == 4)                         % DE/best/2
           FM_ui = FM_bm + (F_weight+F_weight*rand(1,1))*(FM_pm1 + FM_pm2 - FM_pm3 - FM_pm4);    % differential variation
           FM_origin = FM_pm3;
  elseif (I_strategy == 5)                          % DE/rand/2
           FM_ui = FM_pm1+(F_weight+F_weight*rand(1,1))*(FM_pm2-FM_pm3+FM_pm4-FM_pm5);         % differential variation
           FM_origin = FM_pm3;
           FM_ui = FM_popold.*FM_mpo + FM_ui.*FM_mui;     % crossover
  elseif (I_strategy == 6)% DE/target-to-best/p
    for k=1:P
        I_index=floor(100*(rand(1,1))+1);
        I_Indexarray(k)=I_index;                               
        for i=1:I_NP
            if (I_index==i)
               if(k==1)
                  I_indexbest=I_index;
                  m=benchmark_func(FM_popold(I_index,:),code);
               elseif (benchmark_func(FM_popold(I_index,:),code)<m)
                       I_indexbest=I_index;
               end
            end
        end
    end
    for k=1:I_NP                              % population filled with the best member
        FM_bm(k,:) = FM_popold(I_indexbest,:);             % of the last iteration
    end    
    FM_ui = FM_popold + 0.85*(FM_bm-FM_popold) + (F_weight+F_weight*rand(1,1))*(FM_pm1 - FM_pm2);
    FM_ui = FM_popold.*FM_mpo + FM_ui.*FM_mui;
    FM_origin = FM_popold;
    elseif(I_strategy == 7)
       if(rand(1,1)<.30)
       for i=1:I_NP   
        a=[floor(100*rand(1,1)+1) floor(100*rand(1,1)+1) floor(100*rand(1,1)+1)];
        if(i==a(1) | i==a(2) | i==a(3))
        a=[floor(100*rand(1,1)+1) floor(100*rand(1,1)+1) floor(100*rand(1,1)+1)];
        x=abs(benchmark_func(FM_popold(a(1),:),code));
        y=abs(benchmark_func(FM_popold(a(2),:),code));
        z=abs(benchmark_func(FM_popold(a(3),:),code));
        p=x+y+z;
        p1=x/p;
        p2=y/p;
        p3=z/p;
        FM_ui(i,:)=(FM_popold(a(1),:)+FM_popold(a(2),:)+FM_popold(a(3),:))/3+(p2-p1)*(FM_popold(a(1),:)-FM_popold(a(2),:))+(p3-p2)*(FM_popold(a(2),:)-FM_popold(a(3),:))+(p1-p3)*(FM_popold(a(3),:)-FM_popold(a(1),:));
        else
        x=abs(benchmark_func(FM_popold(a(1),:),code));
        y=abs(benchmark_func(FM_popold(a(2),:),code));
        z=abs(benchmark_func(FM_popold(a(3),:),code));
        p=x+y+z;
        p1=x/p;
        p2=y/p;
        p3=z/p;
        FM_ui(i,:)=(FM_popold(a(1),:)+FM_popold(a(2),:)+FM_popold(a(3),:))/3+(p2-p1)*(FM_popold(a(1),:)-FM_popold(a(2),:))+(p3-p2)*(FM_popold(a(2),:)-FM_popold(a(3),:))+(p1-p3)*(FM_popold(a(3),:)-FM_popold(a(1),:));
        end
       end
        FM_ui = FM_popold.*FM_mpo + FM_ui.*FM_mui;
       else
       FM_ui = FM_pm3 + (F_weight+F_weight*rand(1,1))*(FM_pm1 - FM_pm2);   % differential variation
       FM_ui = FM_popold.*FM_mpo + FM_ui.*FM_mui;     % crossover
       end
       else                   
     % either-or-algorithm
      fprintf(1,'Please enter a correct choice');     
  end
  
%-----Optional parent+child selection-----------------------------------------
  
%-----Select which vectors are allowed to enter the new population------------
  for k=1:I_NP
  
     F_tempval= benchmark_func(FM_ui(k,:),code);   % check cost of competitor
      I_nfeval  = I_nfeval + 1;
      if (left_win(F_tempval,F_val(k)) == 1)   
         FM_pop(k,:) = FM_ui(k,:);                    % replace old vector with new one (for new iteration)
         F_val(k)   = F_tempval;                      % save value in "cost array"
      
         %----we update F_bestval only in case of success to save time-----------
         if (left_win(F_tempval,F_bestval) == 1)   
            F_bestval = F_tempval;                    % new best value
            FVr_bestmem = FM_ui(k,:);                 % new best parameter vector ever
         end
      end
   end % for k = 1:NP

  FVr_bestmemit = FVr_bestmem;       % freeze the best member of this iteration for the coming 
                                     % iteration. This is needed for some of the strategies.

%----Output section----------------------------------------------------------


  if (I_refresh > 0)
     if ((rem(I_iter,I_refresh) == 0) | I_iter == 1)
       fprintf(1,'Iteration: %d,  Best: %f,  F_weight: %f,  F_CR: %f,  I_NP: %d   d: %d\n',I_iter,F_bestval,F_weight,F_CR,I_NP,d);
       %var(FM_pop)
       format long e;
       %for n=1:I_D
        %  fprintf(1,'best(%d) = %g\n ',n,FVr_bestmem(n));
       %end
       %I_Indexarray
     end
  end
  

  I_iter = I_iter + 1;
end %---end while ((I_iter < I_itermax) ...
fprintf(1,'FES= %d \n',I_nfeval);
