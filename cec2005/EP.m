clear all;
N=100;
dim=10; %no of dimensions
run=1;%no of runs
num=1; 
%max_fe=100000;
%cut_off=9*10^(-5);
%penalty=20;
a1=sqrt(N);
b1=sqrt(2*a1);
t1=1/b1;
t2=1/sqrt(2*N);
range=repmat([-100 100],dim,1);    % range
iter=1000; %no of iterations
dwnsrch=range(:,1);
upsrch=range(:,2);
FEs = 0;
b=1;
y=0;
q=10;

t=zeros(1,run);
flag=zeros(1,run);

for b=1:run
    
  pop=zeros(N,dim);
  sd_pop=zeros(N,dim)+3;
  sd-off=zeros(N,dim);
  for i=1:N
    for j=1:dim
      pop(i,j)=dwnsrch(j)+rand*(upsrch(j)-dwnsrch(j));
    end
  end
  
  for v=1:iter
  
   rand('state',sum(100*clock));
   off=zeros(N,dim);
    
    for i=1:N 
     a=randn;
     for j=1:dim
       sd_off(i,j)=sd_pop(i,j)*exp(t1*a+t2*randn);
       off(i,j)=pop(i,j)+sd_off(i,j)*randn;
       
           if(off(i,j)>upsrch(j))
             off(i,j)=upsrch(j)-rand;% z(i,j)=upsrch(j)-rand*(upsrch(j)-dwnsrch(j));
           end
           if(off(i,j)<dwnsrch(j))
              off(i,j)=dwnsrch(j)+rand; %z(i,j)=dwnsrch(j)+rand*(upsrch(j)-dwnsrch(j));
           end
     end
    end
         
     for j=1:N
      fitness_pop(j,1)=benchmark_func(pop(j,:),num);
      fitness_pop(j+N,1)=benchmark_func(off(j,:),num);
      FEs=FEs+2;
     end
     
     for k=1:2*N
         win(k,1)=0;
     end
     
     for g=1:N
         for p=1:dim
             comb_pop(g,p)=pop(g,p);
         end
     end
     
     for g=N:2*N
         for p=1:dim
             comb_pop(g,p)=off(g,p);
         end
     end
     
     
     j=0;
     for k=1:2*N
         for d=1:q
             while(k==j)
              j=ceil(2*N*rand);
             end
             if(fitness_pop(k)<fitness_pop(j))
                 win(k,1)=win(k,1)+1;
             end
         end
     end
     
     [win_sort, I]=sort(win,'ascend');
     
          for i=1:N
              for j=1:dim
                  pop(i,j)=comb_pop(I(i),j)
              end
          end
          
          s=sprintf('best fitness in iteration no. %d is = %e',v,y+450);
          disp(s)
  end
end

              
              
       
       
       
       
       
       
       
       
       
       
       