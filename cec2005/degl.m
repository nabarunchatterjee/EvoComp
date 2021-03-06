%DEGL MATLAB Source Code
%Dr. Swagatam Das
%search space defination must be functn specific
clc
clear all
v=100; %no of population vectors
dim=30; %no of dimensions
k=10;   % radius,shud be less (v-1)/2
run=1;%no of runs
num=1; 
fbias=-450;
%max_fe=100000;
%cut_off=9*10^(-5);
%penalty=20;
range=repmat([-100 100],dim,1);    % range
iter=1000; %no of iterations
F=0.8;    % parameter
w=0.45;   % parameter
cr=0.3;   % crossover factor
dwnsrch=range(:,1);
upsrch=range(:,2);
FEs = 0;
b=1;

t=zeros(1,run);
flag=zeros(1,run);

best_vctr=zeros(run,dim);

tmp=0;
pataka=zeros(1,run);
for b=1:run
    
  pop=zeros(v,dim);
  for i=1:v
    for j=1:dim
      pop(i,j)=dwnsrch(j)+rand*(upsrch(j)-dwnsrch(j));
    end
  end
  

   %if(flag(b)>max_fe) % 
    %  break
   %end
    c=zeros(v,1);
   for j=1:v
     c(j,1)=benchmark_func(pop(j,:),num);
     FEs=FEs+1;
%      flag(b)=flag(b)+1;
   end
   
  for a=1:iter
   %k=4+ceil(rand*6);
   rand('state',sum(100*clock));
   [g2,h2]=min(c);
   z1=zeros(v,dim);
   z2=zeros(v,dim);
   z=zeros(v,dim);
   k1=2*k+1;
   d=zeros(k1,1);
   % de mutation AND crossover step 
   for i=1:v
       
       q1=0;
        r1=0;
     while(q1==0||r1==0||q1==r1)
       q1=ceil(k*rand);
       r1=ceil(k*rand);
    end
       if (rand>0.5)
           q1=-q1;
       end
       if (rand>0.5)
           r1=-r1;
       end
       k2=2*k;
       for m=0:k2
         tr=i-k+m;
         tr=circadjst(tr,v);
         d(m+1,1)=c(tr,1);
       end
       [g1,h1]=min(d);
       m1=i-k-1+h1;
       m1=circadjst(m1,v);
       m2=i-q1;
       m2=circadjst(m2,v);
       m3=i-r1;
       m3=circadjst(m3,v);
       for j=1:dim
        %   F=sin(2*pi*rand);
       z1(i,j)=pop(i,j)+F*(pop(m1,j)-pop(i,j))+F*(pop(m2,j)-pop(m3,j));    
       end
       q2=0;
        r2=0;
     while(q2==0||r2==0||q2==r2)
       q2=ceil(v*rand);
       r2=ceil(v*rand);
       %r3=ceil(v*rand);
     end
    for j=1:dim
      %  F=sin(2*pi*rand);
     z2(i,j)=pop(h2,j)+F*(pop(q2,j)-pop(r2,j))+F*(pop(h2,j)-pop(i,j));
    % z2(i,j)=pop(r3,j)+(0.1+0.9*rand)*(pop(q2,j)-pop(r2,j));
    end
     
      b1=ceil(rand*dim);
      %w=1;      
      w=0.5;
      
      for j=1:dim
           if ((sin(2*pi*rand)<=cr)||(j==b1))            % crossover incorporated
           z(i,j)=w*z2(i,j)+(1-w)*z1(i,j);
           if(z(i,j)>upsrch(j))
             %z(i,j)=upsrch(j)-rand;
             z(i,j)=upsrch(j)-rand*(upsrch(j)-dwnsrch(j));
           end
           if(z(i,j)<dwnsrch(j))
             % z(i,j)=dwnsrch(j)+rand; 
             z(i,j)=dwnsrch(j)+rand*(upsrch(j)-dwnsrch(j));
           end
         else
             z(i,j)=pop(i,j);
         end
       end
  
   x=benchmark_func(z(i,:),num);
   FEs=FEs+1;
%    flag(b)=flag(b)+1;
   if(x<=c(i,1))
       pop(i,:)=z(i,:);
       c(i,1)=x;
   end
    end
  if(a==1)
    y=g2;
   end
   if(y>g2)
    y=g2;
   end
  % if(y<cut_off)%stopping condition to check speed
   %   pataka(b)=1;
    %   break
   %end
   if(rem(a,1)==0)
       s=sprintf('best fitness in iteration no. %d is = %f',a,y-fbias);
       disp(s)
   end
   
 end
 
 
 
  x=sprintf('for run no. %d,the best value of fitness is:   ',b);
  disp(x)
  c=sprintf('%e  ',y-360);
  disp(c)
  t(b)=y;
  
  best_vctr(b,:)=pop(h2,:);
  
 end


mean=sum(t)/run;
disp('the average value of error is :  ')
show=sprintf('%e ',mean-360);
disp(show)

sd=std(t-360);
disp('the standard deviation of fitness is :  ')
e=sprintf('%f ',sd);
disp(e)


[o1,o2]=min(t);
best_vctr(o2,:)

disp('the No. of FEs is :  ')
FEs
