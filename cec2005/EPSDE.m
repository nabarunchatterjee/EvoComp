%Differential Evolution
%D--Dimensions
%NP--Population of agents
%Array Cr--Crossover Rate
%Arrays F,K--Scaling Factors
%num--Select number according to function
%bound--assign bound of test function according to function
%Array MS--DE startegy used
%Arrays good_MS, good_F, good_K, good_CR--Database of successful parameters
%fs--No. of function evaluations
%G--Generation
clear all;
%max_runs=10; error_list=zeros(max_runs);
%for runs=1:max_runs
    
D=30; NP=50; 
Cr=0.1:0.1:0.9; KK=0.4:0.1:0.9; FF=0.4:0.1:0.9; %Entire pool of parameters
num=9; bound=5; f_bias=-330;
fs=0;  G=0; MS=zeros(1,NP); CR=zeros(1,NP); F=zeros(1,NP); K=zeros(1,NP); funcval=zeros(1,NP); y=zeros(1,NP); 


%Initialization
X=-bound+rand(NP,D)*2*bound;
X_best=-bound+rand(1,D)*2*bound;
for i=1:NP
   MS(i)=randi(3);
   CR(i)=Cr(randi(9));
   F(i)=FF(randi(6));
   K(i)=KK(randi(6));
end
while fs<=300000  %Checking no. of function evaluations
    
   for i=1:NP
       
        %Mutation
        T1=randi(NP);
        while T1==i
            T1=randi(NP);
        end
        T2=randi(NP);
        while(T2==T1)||(T2==i)
            T2=randi(NP);
        end
        T3=randi(NP);
        while (T3==T2)||(T3==T1)||(T3==i)
            T3=randi(NP);
        end
        T4=randi(NP);
        while(T4==T1)||(T4==T2)||(T4==T3)||(T4==i)
            T4=randi(NP);
        end
        T5=randi(NP);
        while(T5==T1)||(T5==T2)||(T5==T3)||(T5==T4)||(T5==i)
            T5=randi(NP);
        end
        switch MS(i)
            case 1  %DE/rand-to-best/2
                V(i,:)=X(i,:)+K(i)*(X_best-X(i,:))+F(i)*(X(T1,:)-X(T2,:)+X(T3,:)-X(T4,:));
            case 2  %DE/rand/1
                V(i,:)=X(T1,:)+F(i)*(X(T2,:)-X(T3,:));
        end        
        if(MS(i)==3) %DE/current-to-rand/1
           U(i,:)=X(i,:)+1*(X(T1,:)-X(i,:))+F(i)*(X(T2,:)-X(T3,:));
           for j=1:D
	    %Limiting parameters within bound
            if(U(i,j)<-bound)
                U(i,j)=-bound;
            end
            if (U(i,j)>bound) 
                U(i,j)=bound;
            end
           end
        else
        jrand=randi(D);
        for j=1:D
	    %Limiting parameters within bound
            if(V(i,j)<-bound)
                V(i,j)=-bound;
            end
            if (V(i,j)>bound) 
                V(i,j)=bound;
            end
            %Binomial Crossover
            if(rand()<=CR(i)||j==jrand)
                U(i,j)=V(i,j);
            else
                U(i,j)=X(i,j);
            end
        end
       end
       funcval(i)=benchmark_func(X(i,:),num);
       %Selection	
       if(benchmark_func(U(i,:),num)<=funcval(i))
          X(i,:)=U(i,:);
          if(y(i)==0||good_MS(i,y(i))~=MS(i)||good_CR(i,y(i))~=CR(i)||good_F(i,y(i))~=F(i)||good_K(i,y(i))~=K(i))
              y(i)=y(i)+1;
          end
          %Storing successful mutation strategies and parameters
          good_MS(i,y(i))=MS(i);
          good_CR(i,y(i))=CR(i);
          good_F(i,y(i))=F(i);
          good_K(i,y(i))=K(i);
          
       else
           if(rand()<=0.5 && y(i)>=1)  %Asigning strategies and parameters from stored successful pool
               y1=randi(y(i));
               MS(i)=good_MS(i,y1);
               CR(i)=good_CR(i,y1);
               F(i)=good_F(i,y1);
               K(i)=good_K(i,y1);
           else     %Assigning new random strategies
               ms=MS(i); 
               while(MS(i)==ms)
                   MS(i)=randi(3);
               end
               cr=CR(i); 
               while(CR(i)==cr) 
                   CR(i)=Cr(randi(9)); 
               end
               ff=F(i);
               while(F(i)==ff)
                   F(i)=FF(randi(6));
               end
               kk=K(i);
               while(K(i)==kk)
                   K(i)=KK(randi(6));
               end
           end
       end
       fs=fs+2;
   end
    fmin=min(funcval);
    m=find(funcval==fmin,1);
    X_best=X(m,:);
    disp(['Error after ',num2str(G),' generations = ',num2str(fmin-f_bias)]);
    G=G+1;
    
end
%error_list(runs)=fmin-f_bias;
%end
%avg=mean(error_list);
%stand_dev=std(error_list);

        
        
