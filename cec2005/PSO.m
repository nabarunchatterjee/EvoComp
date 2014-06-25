%Script for a basic PSO Algorithm
%For Sphere function     (D=10) set num = 1,
% bound = 100  , dim = 10, Optimum = [ 0 0 0 ... 0 ].
%For Rosenbrock function (D=10) set num = 2,
% bound = 50   , dim = 10, Optimum = [ 1 1 1 ... 1 ].
%For Rastrigin function  (D=10) set num = 3,
% bound = 5.12 , dim = 10, Optimum = [ 0 0 0 ... 0 ].
%For Griewank function   (D=10) set num = 4,
% bound = 600 , dim = 10, Optimum = [ 0 0 0 ... 0 ].
%--------------------------------
max_iter = 150;              %Maximum no. of iterations
dim = 10;
agents = 60;
omega = 0.4;
phi1 = 2.0;
phi2 = 2.0;
x = [];                      %Empty matrix for the agents' positions
v = [];                      %Empty matrix for the agents' velocities
funcval = [];                %Empty matrix for holding the value of the 
num = 5;                     %function at the agents' local best positions
switch num
    case 1
        f = @(x) sum(x.^2);     %Sphere function
        optimum = zeros(1,dim);
        bound = 100;
    case 2
        f = @(x) sum(100*(x(2:dim)-x(1:dim-1).^2).^2 + (1-x(1:dim-1)).^2);
        optimum = ones(1,dim);  %Rosenbrock function
        bound = 50;
    case 3
        f = @(x) sum(x.^2 - 10*cos(2*pi*x) + 10);
        optimum = zeros(1,dim); %Rastrigin function
        bound = 5.12;
    case 4
        indices = 1:dim;
        f = @(x) (1/4000)*sum(x.^2) - prod(cos (x./sqrt(indices))) + 1;
        optimum= zeros(1,dim);
        bound = 600;
        %Griewank function
    case 5
       load(ackley_M_D30.mat);
end
for a = 1:agents
    x(a, : ) = bound * (2*rand(1,dim) - ones(1,dim));%Generating agents
    v(a, : ) = bound * (2*rand(1,dim) - ones(1,dim));%Generating velocities
    funcval(a, : ) = f(x(a, : ));
end
for a = 1:agents
    if min(funcval) == f(x(a, : ))      %Getting the global best postion
        Pgb = x(a, : );                 %Before the start of the algorithm
    end
end
Plb = x;
for k = 1:max_iter
    for m = 1:agents
        v(m, : ) = omega*v(m, : ) + phi1*rand*(Plb(m)-x(m, : )) + phi2*rand*(Pgb - x(m, : ));
        x(m, : ) = x(m, : ) + v(m, : );  %The algorithm
        if f(x(m, : )) <= funcval(m, :)  %Updating local best position
            funcval(m, :) = f(x(m, : )); %if new Plb found   
            Plb(m, : ) = x(m, : );
        end
    end
    if min(funcval)<f(Pgb)
        for a = 1:agents
            if min(funcval) == f(x(a, : ))  %Updating global best position
            Pgb = x(a, : );                 %if new Pgb found
            end
        end
    end
    err = Pgb - optimum;                
    error = sqrt(sum(err.^2));          
    disp(['Error after ',num2str(k),' iterations = ',num2str(error)])
    %keyboard
end