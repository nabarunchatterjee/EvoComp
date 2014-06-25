% Evolutionary Programming using Gaussian Mutation
clear all
%path(path,'D:\Rammohan\EP\DAS\SADE_TEC');

global initial_flag
initial_flag=0;
I_fno=13;
NP=100;
D=30;
fbias=-130;
bound=100;
bound1=3;
bound2=1;

FVr_minbound = -bound1*ones(1,D); 
      FVr_maxbound = bound2*ones(1,D); 
      Lbound       = FVr_minbound;
      Ubound       = FVr_maxbound;


 q=20;
 
Max_Gen=3000;

for runs=1:1
     
  fprintf('Run %d',runs);
  nfeval=0;
  initial_flag=0;


  gen=1;

  pop1=repmat(FVr_minbound,NP,1)+repmat((FVr_maxbound-FVr_minbound),NP,1).*rand(NP,D);
  for i=1:NP
      val1(i,1)=benchmark_func(pop1(i,:),I_fno);
      %val1(i)=feval('DE_benchmark_func',pop1(i),I_fno);
  nfeval=nfeval+1;
  end
  neta1=[(0.8/sqrt(D))*rand(NP,D)].*repmat((FVr_maxbound-FVr_minbound),NP,1);

  best=min(val1);
  NN=0;
    
while(gen<Max_Gen)
     
       sneta=[];
       
       %x=100
      newneta1 = neta1.*[2*rand(NP,D)-1]+neta1; 
      rama=randint(NN,2,[1,NP]);
      new1(1:NN,:)=repmat(pop1(1,:),NN,1)+0.85*(pop1(rama(1:NN,1),:)-pop1(rama(1:NN,2),:));
      new1((NN+1:NP),:)=pop1(NN+1:NP,:)+newneta1(NN+1:NP,:).*cauchy((NP-NN),D);
         
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
    
      s=sprintf('best fitness in iteration no. %d is = %e',gen,bestit-fbias);
   disp(s)
      
end
%eval(['save /staff2/mall/snetacauchy10D/F' int2str(I_fno) '/W_' int2str(runs) ' W;']);
 


Re(runs,:)=[best,nfeval];
end
%eval(['save /staff2/mall/snetacauchy10D/F' int2str(I_fno) '/FUN' int2str(I_fno)]);