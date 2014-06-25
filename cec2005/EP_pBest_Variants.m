clear all;
N=100;
dim=30; %no of dimensions
D=dim;
NP=N;
run=1;%no of runs
num=6;
bound=100;
fbias=390;
%max_fe=100000;
%cut_off=9*10^(-5);
%penalty=20;
a1=sqrt(N);
b1=sqrt(2*a1);
t1=1/b1;
t2=1/sqrt(2*N);
range=repmat([-bound bound],dim,1);    % range
iter=3000; %no of iterations
dwnsrch=range(:,1);
upsrch=range(:,2);
FEs = 0;
b=1;
y=0;
q=10;
si=0.8/sqrt(N);

t=zeros(1,run);
flag=zeros(1,run);

for b=1:run
    
  pop=zeros(N,dim);
  sd_pop=zeros(N,dim)+3;
  
  for i=1:N
    for j=1:dim
      pop(i,j)=dwnsrch(j)+rand*(upsrch(j)-dwnsrch(j));
      %sd_pop(i,j)=si*rand*(upsrch(j)-dwnsrch(j));
    end
  end
  
  for j=1:N
        fitness_pop(j,1)=benchmark_func(pop(j,:),num);
        fitness_pop_only(j,1)=fitness_pop(j,1);
        FEs=FEs+1;
  end    
        
 
  for v=1:iter
  
  
   s=sprintf('best fitness in iteration no. %d is = %f',v,min(fitness_pop_only)-fbias);
   disp(s)
          
   rand('state',sum(100*clock));
   off=zeros(N,dim);
   sd_off=zeros(N,dim);
    
   [s1,f]=sort(fitness_pop_only, 'ascend');
   
   
   
   
      Z=cauchy(NP,D);
      B1=randn(NP,D);
       p1=20;
       
   for i=1:N 
     a=randn;
     R=randn(D,D);
       [Q,S1]=qr(R);
   
     %[s1,f]=sort(fitness_pop_only, 'ascend');
         
     for j=1:dim
       ind=f(ceil(rand*p1));
      
       if(rand<=0.5)
        sd_off(i,j)=sd_pop(ind,j)*exp(t1*a+t2*randn);
        
        off(i,j)=pop(ind,j)+sd_off(ind,j)*randn;      %cauchyrnd(0,1,1);
     
       else  
           sd_off(i,j)=sd_pop(i,j)*exp(t1*a+t2*randn);
            off(i,j)=pop(i,j)+sd_off(i,j)*Z(i,j);  
           %off(i,:)=(R*(pop(i,:)+sd_off(i,:)*Z(1,1)).').';
           
           %off(i,j)=pop(i,j)+sd_off(i,j)*randn;         %cauchyrnd(0,1,1);
       end
       
           if(off(i,j)>upsrch(j))
             off(i,j)=upsrch(j)-rand;
            % z(i,j)=upsrch(j)-rand*(upsrch(j)-dwnsrch(j));
           end
           if(off(i,j)<dwnsrch(j))
              off(i,j)=dwnsrch(j)+rand; %
           %   z(i,j)=dwnsrch(j)+rand*(upsrch(j)-dwnsrch(j));
           end
     end
    end
         
     for j=1:N
      %fitness_pop(j,1)=benchmark_func(pop(j,:),num);
      fitness_pop(j+N,1)=benchmark_func(off(j,:),num);
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
              fitness_pop_only(i,1)=fitness_pop(I(i),1);
          end
          
        for j=1:N
        fitness_pop(j,1)=fitness_pop_only(j,1);
        end    
               
          
  end
  %FEs
  
  x=sprintf('for run no. %d,the best value of error is:   ',b);
  disp(x)
  c=sprintf('%e  ',min(fitness_pop_only)-fbias);
  disp(c)
  t(b)=min(fitness_pop_only)-fbias;
  
  %best_vctr(b,:)=pop(h2,:);
  
 end


mean=sum(t)/run;
disp('the average value of error is :  ')
show=sprintf('%e ',mean);
disp(show)

sd1=std(t-fbias);
disp('the standard deviation of fitness is :  ')
e=sprintf('%f ',sd1);
disp(e)


%[o1,o2]=min(t);
%best_vctr(o2,:)

disp('the No. of FEs is :  ')
FEs

  
   
              
       
       
       
       
       
       
       
       
       
       
       