function f=my_benchmark_func(x,D)

sum=0;
sum1=0;

for i=1:D
    %sum=sum+100*(x(i)^2-x(i+1))^2+(x(i)-1)^2;
    sum=sum+(x(i))^2;
    sum1=sum1+x(i);
end

    f=-20*exp(-0.2*sqrt(0.5*(sum)))-exp(0.5*(cos(2*3.141592*sum1)))+20+exp(1);
%end

%f=sum;