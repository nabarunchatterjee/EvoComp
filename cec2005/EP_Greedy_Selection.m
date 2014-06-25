clear all;
N=100;
dim=10; %no of dimensions
run=1;%no of runs
num=9; 
bound=50;
fbias=-330;
%max_fe=100000;
%cut_off=9*10^(-5);
%penalty=20;
a1=sqrt(N);
b1=sqrt(2*a1);
t1=1/b1;
t2=1/sqrt(2*N);
range=repmat([-bound bound],dim,1);    % range
iter=500; %no of iterations
dwnsrch=range(:,1);
upsrch=range(:,2);
FEs = 0;
b=1;
y=0;
q=25;

t=zeros(1,run);
flag=zeros(1,run);

for b=1:run
    
  pop=zeros(N,dim);
  sd_pop=zeros(N,dim)+3;
  
  for i=1:N
    for j=1:dim
      pop(i,j)=dwnsrch(j)+rand*(upsrch(j)-dwnsrch(j));
    end
  end
  
  for v=1:iter
  
      for j=1:N
        fitness_pop(j,1)=benchmark_func(pop(j,:),num);
        FEs=FEs+1;
      end
      
   rand('state',sum(100*clock));
   off=zeros(N,dim);
   sd_off=zeros(N,dim);
    for i=1:N 
     a=randn;
     [s,f]=min(fitness_pop);
     %[s1,f1]=max(fitness_pop);
     for j=1:dim
       sd_off(i,j)=sd_pop(i,j)*exp(t1*a+t2*randn);
      %off(i,j)=pop(i,j)+sd_off(i,j)*randn;
       if(rand<=0.5)      
        off(i,j)=pop(i,j)+sd_off(i,j)*cauchyrnd(0,1,1);
       %elseif (0.5<rand<=0.7)      
       %   off(i,j)=pop(f1,j)+sd_off(i,j)*cauchyrnd(0,1,1);
        else                
           off(i,j)=pop(i,j)+sd_off(i,j)*cauchyrnd(0,1,1);
       end
       
           if(off(i,j)>upsrch(j))
             off(i,j)=upsrch(j)-rand;% z(i,j)=upsrch(j)-rand*(upsrch(j)-dwnsrch(j));
           end
           if(off(i,j)<dwnsrch(j))
              off(i,j)=dwnsrch(j)+rand; %z(i,j)=dwnsrch(j)+rand*(upsrch(j)-dwnsrch(j));
           end
     end
    end
         
     for j=1:N
      %fitness_pop(j,1)=benchmark_func(pop(j,:),num);
      fitness_off(j,1)=benchmark_func(off(j,:),num);
      FEs=FEs+1;
     end
     
     for k=1:N
         for j=1:dim
             
             if(fitness_off(k)<=fitness_pop(k))
                  pop(k,j)=off(k,j);
                  sd_pop(k,j)=sd_off(k,j);
              end
         end
     end
          
          s=sprintf('best fitness in iteration no. %d is = %e',v,min(fitness_pop)-fbias);
          disp(s)
  end
end

              
              
       
       
       
       
       
       
       
       
       
       
       