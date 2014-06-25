clc
clear all
v=100; %no of population vectors
dim=10; %no of dimensions
D=dim;
run=1;%no of runs
num=1; 
%max_fe=100000;
%cut_off=9*10^(-5);
%penalty=20;
bound=10;
range=repmat([-100 100],dim,1);    % range
fbias=-450;
iter=200; %no of iterations
a1=0.8;    % parameter
a2=0.8;   % parameter
cr=0.9;   % crossover factor
dwnsrch=range(:,1);
upsrch=range(:,2);
FEs = 0;
b=1;
y=0;

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
  

   c=zeros(v,1);
   for j=1:v
     c(j,1)=benchmark_func(pop(j,:),num);
     FEs=FEs+1;
%      flag(b)=flag(b)+1;
   end
   
  for a=1:iter
  
   rand('state',sum(100*clock));
   [g2,h2]=min(c);
   [g3,h3]=max(c);
   for i=1:v
       
       
        q2=0;
        r2=0;
     while(q2==0||r2==0||q2==r2||r2==r3||q2==r3)
       q2=ceil(v*rand);
       r2=ceil(v*rand);
       r3=ceil(v*rand);
     end
    
     %R=randn(D,D);
     %[Q,W]=qr(R);
     %if(rand<=0.5) 
    % z2(i,:)=(Q'*(pop(i,:).')).'+0.8*(pop(q2,:)-pop(r2,:));
    % else
         z2(i,:)=pop(r3,:)+0.8*(pop(q2,:)-pop(r2,:));
     %end
     
     
    % for j=1:dim
       %z2(i,j)=pop(h2,j)+(0.1+0.9*rand)*(pop(q2,j)-pop(r2,j));%+(0.1+0.9*rand)*(pop(h2,j)-pop(i,j));
         %z2(i,j)=pop(i,j)+(0.1+0.9*rand)*(pop(q2,j)-pop(r2,j));
        %z2(i,j)=pop(i,j)+(0.1+0.9*rand)*(pop(q2,j)-pop(r2,j));
       %z2(i,j)=pop(h3,j)+(0.1+0.9*rand)*(pop(q2,j)-pop(r2,j))+(0.1+0.9*rand)*(pop(h2,j)-pop(i,j));
   % end
     
      b1=ceil(rand*dim);
      for j=1:dim
           if ((rand<=cr)||(j==b1))            % crossover incorporated
           z(i,j)=z2(i,j); %w*z2(i,j)+(1-w)*z1(i,j);
           if(z(i,j)>upsrch(j))
             %z(i,j)=upsrch(j)-rand;
             %z(i,j)=upsrch(j)-rand*(upsrch(j)-dwnsrch(j));
           end
           if(z(i,j)<dwnsrch(j))
              %z(i,j)=dwnsrch(j)+rand; %
              %z(i,j)=dwnsrch(j)+rand*(upsrch(j)-dwnsrch(j));
           end
         else
             z(i,j)=pop(i,j);
         end
      end
  
      %for i=1:v
          G(i+a*v-v,:)=z(i,:);
   %   end
       
      
   x=benchmark_func(z(i,:),num);
   FEs=FEs+1;
%    flag(b)=flag(b)+1;
  % if(x<=c(i,1))
      % pop(i,:)=z(i,:);
      % c(i,1)=x;
  % end
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
       s=sprintf('best fitness in iteration no. %d is = %e',a,y-fbias);
       disp(s)
   end
   
 end
 
 
 
  x=sprintf('for run no. %d,the best value of fitness is:   ',b);
  disp(x)
  c=sprintf('%e  ',y-fbias);
  disp(c)
  t(b)=y;
  
  best_vctr(b,:)=pop(h2,:);
  
 end


mean=sum(t)/run;
disp('the average value of error is :  ')
show=sprintf('%e ',mean+330);
disp(show)

sd=std(t+450);
disp('the standard deviation of fitness is :  ')
e=sprintf('%f ',sd);
disp(e)


[o1,o2]=min(t);
best_vctr(o2,:)

disp('the No. of FEs is :  ')
FEs

plot(G(:,1),G(:,2),'k.', 'markersize',5)
hold on;
A1=ceil(min(G(:,1)));
B1=ceil(max(G(:,1)));
A2=ceil(min(G(:,2)));
B2=ceil(max(G(:,2)));
X=transpose([-17:0.5:17]);
Y=transpose([-17:0.5:17]);
for i=1:size(X,1)
    for j=1:size(Y,1)
        Z(i,j)=X(j)^2+Y(i)^2;
    end
end
contour(X, Y, Z);
