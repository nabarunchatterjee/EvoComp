clear all;
N=100;
dim=10; %no of dimensions
run=1;%no of runs
num=16;
bound=5;
fbias=120;
a1=sqrt(N);
b1=sqrt(2*a1);
t1=1/b1;
t2=1/sqrt(2*N);
range=repmat([-bound bound],dim,1);    % range
iter=1000; %no of iterations
dwnsrch=range(:,1);
upsrch=range(:,2);
FEs = 0;
b=1;
y=0;
q=10;
si=0.8/sqrt(N);
NP=N;
D=dim;
alpha=0.4;

t=zeros(1,run);
flag=zeros(1,run);

for b=1:run
    
  pop=zeros(N,dim);
  off=zeros(N,dim);
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
  
  
   s=sprintf('best fitness in iteration no. %d is = %e',v,min(fitness_pop_only)-fbias);
   disp(s)
          
   %rand('state',sum(100*clock));
   off=zeros(N,dim);
   sd_off=zeros(N,dim);
    
   ZN=randn(N,dim);
   Z=cauchy(N,dim);
   B=randn(N,dim);
   
   
   
   %p1=20;
   %p1=ceil(0.5*(N-N*((v-1)/iter)));
   
   %for i=1:N 
     [s1,f]=sort(fitness_pop_only, 'ascend');
     %[s1,f1]=min(fitness_pop);
     
     %a=randn;
     
     %x1=ceil(0.5*rand);
     
      Z=cauchy(NP,D);
      B1=randn(NP,D);
         p1=10;
          %p1=ceil(0.4*(NP-NP*((gen-1)/Max_Gen)));
       for i=1:NP
          %R=rand(D,D);
          %[Q,S1]=qr(R);
          a=randn;
          ind=f(ceil(rand*p1));
          sd_off(i,:)=sd_pop(i,:).*exp(t1*a+t2*B(i,:));
          ind=f(ceil(rand*p1));
          %R=rand(D,D);
          %[Q,S1]=qr(R);
          
          if(rand<=0.5)
              %for j=1:D
               off_temp(i,:)=pop(i,:)+sd_off(i,:).*Z(i,:);
               %off(i,:)=alpha*pop(i,:)+(1-alpha)*pop(ind,:);
              %end
         else  
            %for j=1:D
               off(i,:)=pop(i,:)+sd_off(i,:).*Z(i,:);
              %off(i,:)=(Q*pop(i,:).').'+sd_off(i,:).*Z(i,:);
         end
      end
      
                
                
     %for j=1:dim
            
      % if(rand<=0.5)
      %  sd_off(i,j)=sd_pop(i,j)*exp(t1*a+t2*B(i,j));
       % off(i,j)=pop(ind,j)+sd_off(i,j)*Z(i,j);
     
     %  else  
      %     sd_off(i,j)=sd_pop(i,j)*exp(t1*a+t2*B(i,j));
      %     off(i,j)=pop(i,j)+sd_off(i,j)*Z(i,j);
    %   end
       for j=1:dim
           if(off(i,j)>upsrch(j))
             off(i,j)=upsrch(j);
            
           end
           if(off(i,j)<dwnsrch(j))
              off(i,j)=dwnsrch(j); %
         
           end
     end
   
         
     for j=1:N
     
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
e=sprintf('%e ',sd1);
disp(e)


%[o1,o2]=min(t);
%best_vctr(o2,:)

disp('the No. of FEs is :  ')
FEs

  
   
              
       
       
       
       
       
       
       
       
       
       
       