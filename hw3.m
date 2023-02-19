%% HW3 Team Homework SI - Modeling
%% AMSC420 - Jiatong Liang and Adam Levav
%% Question 1 Preprocessing Data

mytable = readtable('project8_data.xlsx');
x = mytable{2:3,13:end};
x = transpose(x);

%% Setting up the data
t0 = -1;
Tmax = 119;

for i = 1:size(x, 1)
    if x(i,1) >= 5
        t0 = i;
        break
    end
end
t0
I = zeros(120,1);

% Here we are copying over 119 ADDITIONAL days so there's 120 days in total
for j = 0:119
    I(j+1, 1) = x(t0 + j, 1);
end

% Max population is pulled from Population column
Nmax = 236842

% Nmin is maximum infected population based on our data of 120 entries
Nmin = 1 + I(120, 1)

%% Problem 1 part a)
N = Nmin;
coeff = 6/(Tmax*(Tmax+1)*(2*Tmax + 1));
summation = 0;

% Because of MATLAB indexing, we want to start with I(2) and end with
% I(120) for the summation since the lecture has indexing with I(0) but in
% MATLAB that's equivalent to I(1)
for t = 1:Tmax
    summation = summation + t*log((I(t+1)/I(1)) * ((N-I(1))/(N-I(t+1))));
end
% This calculates betahat
betahat = coeff * summation
J = 0;

% This will calculate J(beta, N)
for t = 1:(Tmax+1)
    J = J + abs(betahat*t - log(I(t)/(N-I(t))) + log(I(1)/(N-I(1))))^2;
end

J
I_fctn = @(x) (N*I(1)) / (I(1) + (N-I(1))*exp(-1*betahat*x));
fplot(I_fctn,[0 120]);

%% Problem 1 part a)
N = Nmax;
coeff = 6/(Tmax*(Tmax+1)*(2*Tmax + 1));
summation = 0;

% Because of MATLAB indexing, we want to start with I(2) and end with
% I(120) for the summation since the lecture has indexing with I(0) but in
% MATLAB that's equivalent to I(1)
for t = 1:Tmax
    summation = summation + t*log((I(t+1)/I(1)) * ((N-I(1))/(N-I(t+1))));
end
% This calculates betahat
betahat = coeff * summation
J = 0;

% This will calculate J(beta, N)
for t = 1:(Tmax+1)
    J = J + abs(betahat*t - log(I(t)/(N-I(t))) + log(I(1)/(N-I(1))))^2;
end

J
I_fctn = @(x) (N*I(1)) / (I(1) + (N-I(1))*exp(-1*betahat*x));
fplot(I_fctn,[0 120]);

%% Problem 1 part b)
J = zeros(1,1);
Jold = inf;
N = 1 + I(120);
a = 6/(Tmax*(Tmax+1)*(2*Tmax + 1));
summation1 = 0;
summation2 = 0;
    
for t = 1:Tmax
    summation2 = summation2 + t*log((I(t+1)/I(1)) * ((N-I(1))/(N-I(t+1))));
end
    
for t = 1:Tmax
    summation1 = summation1 + abs(log((I(t+1)/I(1)) * ((N-I(1))/(N-I(t+1)))))^2;
end

JN = summation1 - a*(summation2^2);
J(1,1) = JN;
i = 2;

while JN < Jold
    Jold = JN;
    N = N+1
    summation1 = 0;
    summation2 = 0;
    
    for t = 1:Tmax
        summation2 = summation2 + t*log((I(t+1)/I(1)) * ((N-I(1))/(N-I(t+1))));
    end
    
    for t = 1:Tmax
        summation1 = summation1 + abs(log((I(t+1)/I(1)) * ((N-I(1))/(N-I(t+1)))))^2;
    end

    JN = summation1 - a*(summation2^2);
    J(i,1) = JN;
    i = i + 1;
end

summation = 0;
for t = 1:Tmax
    summation = summation + t*log((I(t+1)/I(1)) * ((N-I(1))/(N-I(t+1))));
end

plot(J);
N
betahat = a*summation

JBN = 0;

% This will calculate J(beta, N)
for t = 1:(Tmax+1)
    JBN = JBN + abs(betahat*t - log(I(t)/(N-I(t))) + log(I(1)/(N-I(1))))^2;
end

JBN

I_fctn = @(x) (N*I(1)) / (I(1) + (N-I(1))*exp(-1*betahat*x));
fplot(I_fctn,[0, 120]);

%% Problem 1 part b) part iv
J = zeros(1,1);
i = 1;
N = Nmin;
JN = 0;

while N < Nmax
    summation1 = 0;
    summation2 = 0;
    
    for t = 1:Tmax
        summation2 = summation2 + t*log((I(t+1)/I(1)) * ((N-I(1))/(N-I(t+1))));
    end
    
    for t = 1:Tmax
        summation1 = summation1 + abs(log((I(t+1)/I(1)) * ((N-I(1))/(N-I(t+1)))))^2;
    end

    JN = summation1 - a*(summation2^2);
    J(i,1) = JN;
    i = i + 1;
    N = N + 1;
end

plot(J)
min(J)

%% Problem 1 part c
new_I = zeros(1,1);
i = 1;
N = Nmin;
coeff = 6/(Tmax*(Tmax+1)*(2*Tmax + 1));

while N < Nmax
    summation = 0;
    X(i) = N;
    % Because of MATLAB indexing, we want to start with I(2) and end with
    % I(120) for the summation since the lecture has indexing with I(0) but in
    % MATLAB that's equivalent to I(1)
    for t = 1:Tmax
        summation = summation + t*log((I(t+1)/I(1)) * ((N-I(1))/(N-I(t+1))));
    end
    % This calculates betahat
    betahat = coeff * summation;
    
    summation = 0;
    for t = 0:Tmax
        summation = summation + abs(I(t+1) - (N*I(1))/(I(1) + (N-I(1))*exp(-1*betahat*t)))^2;
    end

    new_I(i,1) = summation;
    i = i + 1;
    N = N + 1;
end

plot(X, new_I)
min(new_I)

disp(2)