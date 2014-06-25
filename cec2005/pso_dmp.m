%% Particle Swarm Optimization with DIfference Mean Based Optimization
function[time]=pso_dmp(~,~)
clc;
clear all
close all
%% Input Requirements
%% Functions
% 1. Sphere
% 2. Schwefels
% 3. Rosenbrock
% 4. Noise
% 5. Schwefel_multimodal
% 6. Rastrigin_multimodal
% 7. Ackley_mutimodal
% 8. Griewank_multimodal
% 9. Generalized
% 10. Penalized
%%
fun_num=input('enter func no.:'); %see above for function details
%% Declaration of Variables
d=50; %dimentionality of the problem
np=1; %no. of subpopulations
r=50; %number of indivduals in a subpopulation
pop=zeros(r,d,np); % population initialization
lbest=zeros(np,d); %local best
row_index_min=zeros(1,np); % index to determine local best
fitness=zeros(np,r); % matrix showing the functional value for each individual in a population
stagnate=zeros(np,r); 
lbestfit=NaN(1,np); % since pso dmp is a subpopulation based approach lbest denotes the local best for each subpopulation
iter=1;
maxiter=4000;
xaxis=zeros(1,floor(maxiter/500)+1); 
yaxis=zeros(1,floor(maxiter/500)+1); 
%% Info function returns the values associated for a particular function
[ ub,lb,range] = info( fun_num);
%% Initialization
for i=1:np
        pop(:,:,i)=lb+(ub-lb)*rand(r,d); % lb=lower bound and ub=upper bound
end
%% initializing velocity matrix
v=zeros(r,d,np);
vmax=0.2*(range); % defining maximum velocity
%% initializing fitness/functional value of initialized variables
for i=1:np
    for j=1:r
        [fitness(i,j)]=psobench(pop(j,:,i),d,fun_num);%fitness=[np x r]
    end
end
%% defining lbest
for i=1:np
row_index_min(i)=1;lbest(i,:)=pop(1,:,i);lbestfit(i)=fitness(i,1);
for j=1:r
    if(fitness(i,row_index_min(i))>fitness(i,j))
        row_index_min(i)=j;%defining index for local best for each subpopulations
        lbest(i,:)=pop(j,:,i);%defining lbest
        lbestfit(i)=fitness(i,j);%storing lbest fitness values to determine global best later
    end
end
end
%% defining pbest
pbest=pop; %initially pbest=particle itself
pbestfit=fitness;
%%
tic
while(iter<=maxiter)
 w=(0.9-0.5*(iter/maxiter)); %inertia weight
 %% Acceleration coefficients
 c1=0.5+0.5*(exp(-iter/500))+1.4*(sin(iter)/30); 
 c2=1+1.4*(1-exp(-iter/500))+1.4*(sin(iter)/30);
 %%
    for i=1:np
        for j=1:r
%% updating velocity
   v(j,:,i)=w*v(j,:,i)+c1*rand(1,d).*(pbest(j,:,i)-pop(j,:,i))+c2*rand(1,d).*(lbest(i,:)-pop(j,:,i));
%% limit velocity vector within predefined bounds        
            for k=1:d
            v(j,k,i)=min(vmax,max(-vmax,v(j,k,i)));
            end
 %% particle update        
            pop(j,:,i)=pop(j,:,i)+ v(j,:,i); 
 %% bound each dimension of particles within specified ower and upper bounds          
            for p=1:d
                pop(j,p,i)=min(ub,max(lb,pop(j,p,i)));
            end
 %% updating fitness
            [fitness(i,j)]=psobench(pop(j,:,i),d,fun_num);
 %% updating pbest
            if(pbestfit(i,j)>fitness(i,j))%if fitness deteriorates or improves accordingly pbest changes
                pbest(j,:,i)=pop(j,:,i);
                pbestfit(i,j)=fitness(i,j);
                stagnate(i,j)=0;
            else stagnate(i,j)=stagnate(i,j)+1;
            end
            if(j~=row_index_min(i))
                %% aging
            if(stagnate(i,j)>10) % checking aging threshold
                pop(j,:,i)=lb+(ub-lb)*rand; 
               [fitness(i,j)]=psobench(pop(j,:,i),d,fun_num);
                pbest(j,:,i)=pop(j,:,i);
                pbestfit(i,j)=fitness(i,j);
                stagnate(i,j)=0;
            end
            end
        end
%% Difference Mean based Perturbation
% Declaration of variables for DMP
best=pop(row_index_min(i),:,i);
bestlocmean=mean(best);  %calculate the dimensional mean of the best individual
% DMP
for k=1:r
if(k~=row_index_min(i))
    ind_mean=mean(pop(k,:,i));  %calculate the dimensional mean of the selected individuals
    diff_mean=(bestlocmean-ind_mean);   %calculate difference mean
    rand_vec=rand(1,d);  %declaration of a vector whose dimensional values belong to a standard
                                                                            %  normal distributioni.e. N(0,1)
    norm_rand=norm(rand_vec);   %calculate the distance of the random vector from origin
   unit_rand=(rand_vec/norm_rand);               %calculate the randomly directed unit vector
    pop(k,:,i)=pop(k,:,i)+diff_mean.*unit_rand;        %difference mean based perturbation of the selected individual
    [fitness(i,k)]=psobench(pop(k,:,i),d,fun_num); % Update fitness of the respective individual
    %% updating pbest
            if(pbestfit(i,k)>fitness(i,k))%if fitness deteriorates or improves accordingly pbest changes
                pbest(k,:,i)=pop(k,:,i);
                pbestfit(i,k)=fitness(i,k);
            end
end
end
% End of DMP
    end
%% Updating local best of each sub population
for i=1:np
row_index_min(i)=1;
for j=1:r
    if(pbestfit(i,row_index_min(i))>pbestfit(i,j))
        row_index_min(i)=j;
        lbest(i,:)=pbest(j,:,i);
        lbestfit(i)=pbestfit(i,j);
    end
end
end
%% Updating  Global best among all the sub populations
minfunc=lbestfit(1);index=row_index_min(1);note=1;
for j=1:np
    if(minfunc>lbestfit(i))
        minfunc=lbestfit(i);
        index=row_index_min(j);
        note=j;
    end
end
%% Result Section
Best_obj_func=(minfunc);
if(iter==1)
    xaxis=iter;
    yaxis=Best_obj_func;
end
if(mod(iter,500)==0)
fprintf('\n The Best objective functional value in %dth iteration is:%e',iter,Best_obj_func);
yaxis=[yaxis,Best_obj_func];
xaxis=[xaxis,iter];
end
iter=iter+1;
end
clc; 
 %% Final output
  disp(sprintf('\n\t      -----------------------------------------Final Output--------------------------------------\n\tThe final optimized value of function %d obtained by PSO-DMP \n\tafter %d iterations is      :     %15.8e', fun_num,maxiter, Best_obj_func)); 
 toc
%% Plotting
semilogy(xaxis,yaxis,'-rs','LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',10);
xlabel('Generations')
ylabel('Best objective functional value in log scale')
grid on
title(sprintf('Convergence graph for function %d obtained by PSO-DMP', fun_num));
%%
end

%% Benchmark
function [ f] = psobench( x,d,fun_num)
if nargin<3
    fprintf('NOT ENOUGH ARGUEMENTS');    
end
switch fun_num
    case 1
      f=sphere(x,d);  
    case 2
    f=schwefel(x,d);
    case 3
    f=rosenbrock(x,d);
    case 4
    f=noise(x,d);
    case 5
    f=schwefel_mult(x,d);
    case 6
    f=rastrigin_mult(x,d);
    case 7
    f=ackley_mult(x,d);
    case 8
    f=griewank_mult(x,d);
    case 9
    f=generalized_mult(x,d);
    case 10
    f=penalised_mult(x,d);
    otherwise
        error('Kindly check the function number.');
end
end

%% Functions
%% Sphere
function [ f ] = sphere( x,d)
f=0;
for i=1:d
f=f+(x(i))^2;
end
end
%% Schwefel
function [ f ] = schwefel(x,d)
sum=0;prod=1;
for i=1:d
    sum=sum+abs(x(i));
    prod=prod*abs(x(i));
end
f=sum+prod;
end
%% Rosenbrock
function [ val ] = rosenbrock( x,d )
val=0;
for i=1:(d-1)
val=val+100*(x(i+1)-x(i)^2)^2+(-x(i)+1)^2;
end
end
%% Noise
function [ f ] = noise( x,d )
sum=0;
for i=1:d
    sum=sum+d*x(i)^4;
end
f=sum+rand;
end
%% Schwefel_multimodal
function [ f ] = schwefel_mult( x,d )
sum=0;
for i=1:d
    sum=sum+x(i)*sin(sqrt(abs(x(i))));
end
f=418.9829*d-sum;    
end
%% Rastrigen_multimodal
function [ val] =rastrigin_mult( x,y )
val=0;
for i=1:y
    val=val+(x(i)^2-10*cos(2*pi*x(i))+10);
end
end
%% Ackley_multimodal
function [ f ] = ackley_mult( x,d )
sum1=0;sum2=0;
for i=1:d
    sum1=sum1+x(i)*x(i);
    sum2=sum2+cos(2*pi*x(i));
end
f=-20*exp(-0.2*sqrt(sum1/d))-exp(sum2/d)+20+exp(1);
end
%% Griewank_multimodal
function [ val ] = griewank_mult(x,y)
d=0;prod=1;
for i=1:y
d=d + (x(i)^2);
prod=prod*cos(x(i)/(sqrt(i)));
end
val=(d/4000)-prod+1;
end
%% Generalized_multimodal
function [ f] = generalized_mult( x,d )
sum1=0;sum2=0;u=zeros(1,d);
for i=1:d
y(i)=((x(i)+1)/4)+1;
if(x(i)>10)
u(i)=100*(x(i)-10)^4;
else
    if(x(i)<-10)
        u(i)=100*(-x(i)-10)^4;
    else
        u(i)=0;
    end
end
sum2=sum2+u(i);
end
for i=1:d-1
    sum1=sum1+(y(i)-1)^2*(1+10*(sin(pi*y(i+1)))^2);
end
f=(pi/d)*((10*(sin(pi*y(1))^2)+sum1+(y(d)-1)^2))+sum2;
end
%% Penalized_multimodal
function [ f ] = penalised_mult( x,d )
sum1=0;sum2=0;u=zeros(1,d);
for i=1:d
if(x(i)>5)
u(i)=100*(x(i)-5)^4;
else
    if(x(i)<-5)
        u(i)=100*(-x(i)-5)^4;
    else
        u(i)=0;
    end
end
sum2=sum2+u(i);
end
for i=1:d-1
    sum1=sum1+(x(i)-1)^2*(1+(sin(3*pi*x(i+1)))^2);
end
f=(1/10)*((sin(3*pi*x(1))^2)+sum1+((x(d)-1)^2)*(1+sin(2*pi*x(d))^2))+sum2;
end
%%
%% information about all functions
function [ ub,lb,range] = info( fun_num)
if nargin<1
    fprintf('NOT ENOUGH ARGUEMENTS');    
end
switch fun_num
    case 1
    ub=50;
    lb=-100;
    range=200;
    case 2
   ub=5;
    lb=-10;
    range=20;
    case 3
     ub=10;
    lb=-10;
    range=20;
   	case 4
    ub=0.64;
    lb=-1.28;
    range=1.28*2;
   case 5
     ub=500;
    lb=-500;
    range=1000;
   case 6
    ub=2;
    lb=-5.12;
    range=5.12*2;
   case 7
     ub=16;
    lb=-32;
    range=32*2;
   case 8
     ub=200;
    lb=-600;
    range=1200;
    case 9
    ub=25;
    lb=-50;
    range=100;
   case 10
    ub=25;
    lb=-50;
    range=100;
    otherwise
        error('Kindly check the function number.');
end
end
%% end

