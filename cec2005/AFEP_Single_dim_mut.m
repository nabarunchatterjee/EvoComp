% Evolutionary Programming using Gaussian Mutation
clear all
%path(path,'D:\Rammohan\EP\DAS\SADE_TEC');

global initial_flag
initial_flag=0;
I_fno=6;
NP=100;
D=30;
fbias=390;%-310;
bound=100;
bound1=3;
bound2=1;
q=30;

FVr_minbound = -bound*ones(1,D); 
      FVr_maxbound = bound*ones(1,D); 
      Lbound       = FVr_minbound;
      Ubound       = FVr_maxbound;


 
Max_Gen=1500;

for runs=1:1
     
  fprintf('Run %d',runs);
  nfeval=0;
  initial_flag=0;


  gen=1;

  pop1=repmat(FVr_minbound,NP,1)+repmat((FVr_maxbound-FVr_minbound),NP,1).*rand(NP,D);
  for i=1:NP
      val1(i,1)=benchmark_func(pop1(i,:),I_fno);
      %val1(i,1)=my_benchmark_func(pop1(i,:),D);
      %val1(i)=feval('DE_benchmark_func',pop1(i),I_fno);
  nfeval=nfeval+1;
  end
  neta1=[(0.8/sqrt(D))*rand(NP,D)].*repmat((FVr_maxbound-FVr_minbound),NP,1);

  best=min(val1);
  NN=0;
    
while(gen<Max_Gen)
     
       sneta=[];
       [s,f]=sort(val1, 'ascend');
       %x=100
      
       newneta1 = neta1.*[2*rand(NP,D)-1]+neta1; 
      
      
      
      rama=randint(NN,2,[1,NP]);
      new1(1:NN,:)=repmat(pop1(1,:),NN,1)+(0.4+0.4*rand)*(pop1(rama(1:NN,1),:)-pop1(rama(1:NN,2),:));
      %new1((NN+1:NP),:)=pop1(NN+1:NP,:)+newneta1(NN+1:NP,:).*cauchy((NP-NN),D);
         
      
      
      Z=cauchy(NP,D);
      B1=randn(NP,D);
          p1=10;
          p2=ceil(NP*((gen-1)^0.92/Max_Gen));
      
       %SP=((Max_Gen-gen)/Max_Gen); 
      if(p2>1)
       for i=NN+1:p2
          ind=f(ceil(rand*p1));
          xs=Z(1,1);
          %for(i=NN:30)
          %if (rand<=0.5)
              %for j=1:D
                k11=ceil(rand*D);
                for j=1:D
                   if (j==k11)
                     new1(i,:)=pop1(ind,:)+newneta1(ind,:).*B1(i,:);
                    else 
                       new1(i,j)=pop1(i,j);
                   end
               end
       end
      end
         %else  
            %for j=1:D
               for(i=p2+1:NP) 
                  
                   if(rand<=0.5)
                    k2=ceil(rand*NP);
                    for j=1:D
                    if (j==k2)
                     new1(i,j)=pop1(i,j)+newneta1(i,j).*Z(i,j);
                     else 
                        new1(i,j)=pop1(i,j);
                    end
                 end
                  else 
                      % ind1=f(ceil(rand*1));
                       i1=0;
                       j1=0;
                       k1=0;
                       while(i1==j1|j1==k1|k1==i1)
                       i1=ceil(rand*100);
                       %ind1=f(ceil(rand*10));
                       j1=ceil(rand*100);
                       k1=ceil(rand*100);
                       end
                       new1(i,:)=pop1(i1,:)+(0.6+rand*0.4)*(pop1(j1,:)-pop1(k1,:));
                       
                   end
                  %new1(i,:)=(Q*pop1(i,:).').'+newneta1(i,:).*Z(i,:);
            %end
        
      end
      
         %new1((1:NP),:)=pop1(1:NP,:)+newneta1(1:NP,:).*cauchy(NP,D);
      
      
      for i=1:NP
          for j=1:D
              if(new1(i,j)<Lbound(j))
                  new1(i,j)=Lbound(j);
              end
              if(new1(i,j)>Ubound(j))
                  new1(i,j)=Ubound(j);
              end
          end
      end
      
    for i=1:NP
           newval1(i,1)=benchmark_func(new1(i,:),I_fno);
           %newval1(i,1)=my_benchmark_func(new1(i,:),D);
            nfeval=nfeval+1;
      end
      
      r=rem(gen,10);
      if(r==0)
        r=10;
      end  
      
     Totpop1=[pop1;new1;];
     Totneta1=[neta1;newneta1];
     Totval1=[val1;newval1];
     win1=zeros(1,2*NP);
     
     
       for kk=1:NP*2
           win1(kk)=0;
           for i=1:q
               competitor1=ceil(rand(1)*NP*2);
               win1(kk)=win1(kk)-(~(Totval1(competitor1)<Totval1(kk)));
               
           end
       end

       [sort_chd1,a] = sortrows([win1',Totpop1,Totval1,Totneta1],1);
       a=a(1:NP);
       pop1=sort_chd1(1:NP,2:(D+1));
       val1=sort_chd1(1:NP,(D+2));
       neta11=sort_chd1(1:NP,(D+3):(D*2+2));
       win1=sort_chd1(1:NP,1)';
      
     if(size(find(a>NP),1)~=0  )
       sneta=[sneta;neta11(find( a>NP),:)];
    else
       sneta=[sneta;(0.8/sqrt(D))*rand(NP,D).*repmat((FVr_maxbound-FVr_minbound),NP,1)]; 
     end
    
    eval(['sneta' '_' num2str(r) '=sneta;']);
    
                       
        L=[];
        for ind=1:min(gen,10)
           L=[L;eval(['sneta_' num2str(ind)])];
        end
           L(~any(L,2),:) = [];
           l=mean(L,1);                 
           neta1=repmat(l,NP,1);
        
    bestit=min(val1);
    
    if(bestit<best)
        best=bestit;
    end
   
     W(gen,1)=best;
     W(gen,2)=nfeval;
     
      gen=gen+1;
    
      s=sprintf('best fitness in iteration no. %d is = %e, p2=%d',gen,bestit-fbias,p2);
   disp(s)
      
end
%eval(['save /staff2/mall/snetacauchy10D/F' int2str(I_fno) '/W_' int2str(runs) ' W;']);
 


Re(runs,:)=[best,nfeval];
end
%eval(['save /staff2/mall/snetacauchy10D/F' int2str(I_fno) '/FUN' int2str(I_fno)]);