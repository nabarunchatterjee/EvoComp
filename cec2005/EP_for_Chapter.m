clear all;
N=50;
dim=2; %no of dimensions
run=1;%no of runs
num=22;
bound=10;
bound1=10;
bound2=-10;
fbias=360;
%max_fe=100000;
%cut_off=9*10^(-5);
%penalty=20;
a1=sqrt(N);
b1=sqrt(2*a1);
t1=1/b1;
t2=1/sqrt(2*N);
range=repmat([-bound bound],dim,1);    % range
iter=25; %no of iterations
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
  
  for i=1:N
    for j=1:dim
      pop(i,j)=dwnsrch(j)+rand*(upsrch(j)-dwnsrch(j));
    end
  end
  
  X=[-10:0.5:10];
Y=[-10:0.5:10];
for i=1:size(Y,2)
     for j=1:size(X,2)
      S(i,j)=-20*exp(-0.2*sqrt(0.5*(X(j)^2+Y(i)^2)))-exp(0.5*(cos(2*3.141592*X(j))+cos(2*3.141592*Y(i))))+20+exp(1);  
     end
   end
   contour(X,Y,S);
   hold on
   plot(pop(:,1),pop(:,2),'k.');
   %hold on;
   figure(1);
  for v=1:iter
  
      if(mod(v,5)==0)
      contour(X,Y,S);
      hold on
      plot(pop(:,1),pop(:,2),'k.');
      figure(v);
      end
      for j=1:N
        fitness_pop(j,1)=my_benchmark_func(pop(j,:),dim);
        fitness_pop_only(j,1)=fitness_pop(j,1);
        FEs=FEs+1;
      end
      
   rand('state',sum(100*clock));
   off=zeros(N,dim);
   sd_off=zeros(N,dim);
    for i=1:N 
     a=randn;
     [s,f]=min(fitness_pop_only);
     %[s1,f1]=max(fitness_pop);
     for j=1:dim
       sd_off(i,j)=sd_pop(i,j)*exp(t1*a+t2*randn);
        off(i,j)=pop(i,j)+sd_off(i,j)*randn;%cauchy(1,1);
       
       
       
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
      fitness_pop(j+N,1)=my_benchmark_func(off(j,:),dim);
      FEs=FEs+1;
     end
     
     for k=1:2*N
         win(k,1)=0;
     end
     
     for g=1:N
         for p=1:dim
             comb_pop(g,p)=pop(g,p);
             comb_sd(g,p)=sd_pop(g,p);
         end
     end
     
     for g=N+1:2*N
         for p=1:dim
             comb_pop(g,p)=off(g-N,p);
             comb_sd(g,p)=sd_off(g-N,p);
         end
     end
     
     
     j=1;
     for k=1:2*N
         for d=1:q
             %while(k==j)
             % j=ceil(2*N*rand);
            % end
             if(fitness_pop(k)<fitness_pop(ceil(2*N*rand)))
                 win(k,1)=win(k,1)+1;
             end
         end
     end
     
     [win_sort, I]=sort(win,'descend');
     
          for i=1:N
              for j=1:dim
                  pop(i,j)=comb_pop(I(i),j);
                  sd_pop(i,j)=comb_sd(I(i),j);
              end
          end
          
          s=sprintf('best fitness in iteration no. %d is = %e',v,min(fitness_pop_only));
          disp(s)
  end
end

              
              
       
       
       
       
       
       
       
       
       
       
       

            
