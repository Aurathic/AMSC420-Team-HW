%% HW3 Team Homework SI - Modeling
%% AMSC420 - Jiatong Liang and Adam Levav
%% Question 1 Preprocessing Data

% Reading in the Excel file as a table then converting it into a matrix
mytable = readtable('project8_data.xlsx');
x = mytable{2:3,13:end};
x = transpose(x);

%% Setting up the data
% The purpose of t0 was to simply figure out the first occurrence where the
% number of infected is greater than 5
t0 = -1;
Tmax = 119;

% Iterate through all the rows of our original data
for i = 1:size(x, 1)
    if x(i,1) >= 5
        t0 = i;
        break
    end
end
t0

% vector I will store the 119 days of infection after the day t0
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
% Here we run Algorithm 1 for Nmin
N = Nmin;
% We set up the coefficient for betahat
coeff = 6/(Tmax*(Tmax+1)*(2*Tmax + 1));
summation = 0;

% Because of MATLAB indexing, we want to start with I(2) and end with
% I(120) for the summation since the lecture has indexing with I(0) but in
% MATLAB that's equivalent to I(1)
% We use for loop to calculate the summation from 1 to Tmax
for t = 1:Tmax
    summation = summation + t*log((I(t+1)/I(1)) * ((N-I(1))/(N-I(t+1))));
end
% This calculates betahat
betahat = coeff * summation
J = 0;

% This will calculate J(beta, N)
for t = 0:Tmax
    J = J + abs(betahat*t - log(I(t+1)/(N-I(t+1))) + log(I(1)/(N-I(1))))^2;
end

J
% I_fctn is the function representing our predicted value of the infection
% count based on the SI model 
figure(1)
I_fctn = @(x) (N*I(1)) / (I(1) + (N-I(1))*exp(-1*betahat*x));
fplot(I_fctn,[0 120]);
hold on
% I is the true infection count
plot(I)
hold off

%% Problem 1 part a)
N = Nmax;
% Resetting summation value
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
for t = 0:Tmax
    J = J + abs(betahat*t - log(I(t+1)/(N-I(t+1))) + log(I(1)/(N-I(1))))^2;
end

J
% I_fctn is the function representing our predicted value of the infection
% count based on the SI model 
figure(2)
I_fctn = @(x) (N*I(1)) / (I(1) + (N-I(1))*exp(-1*betahat*x));
fplot(I_fctn,[0 120]);
hold on
% I holds the true infection count
plot(I)
hold off

%% Problem 1 part b)
% We redefine variable J to store the values of J(N)
J = zeros(1,1);
Jold = inf;
N = 1 + I(120);
a = 6/(Tmax*(Tmax+1)*(2*Tmax + 1));

% We need to calculate the first value of J(N) before we run the while loop
% summation1 is the first leftmost summation in the equation
summation1 = 0;
% summation2 is the second rightmost summation in the equation
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
    N = N+1;
    % We have to reset the summations everytime before we calculate
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

figure(3)
plot(J);
N
betahat = a*summation

JBN = 0;
% This will calculate J(beta, N)
for t = 0:Tmax
    JBN = JBN + abs(betahat*t - log(I(t+1)/(N-I(t+1))) + log(I(1)/(N-I(1))))^2;
end
JBN

% I_fctn is the function representing our predicted value of the infection
% count based on the SI model 
figure(4)
I_fctn = @(x) (N*I(1)) / (I(1) + (N-I(1))*exp(-1*betahat*x));
fplot(I_fctn,[0 120]);
hold on
% I holds the true infection count
plot(I)
hold off

%% Problem 1 part b) part iv
J = zeros(1,1);
i = 1;
N = Nmin;
JN = 0;

% Here we run through Nmin to Nmax
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

figure(5)
plot((1:3000), J(1:3000))
min(J)

%% Problem 1 part c
% Here we store values for the ideal objective function
ideal_I = zeros(1,1);
i = 1;
N = Nmin;

% Running from Nmin to Nmax
while N < Nmax
    %% First we need to calculate betahat for each N
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
    
    % Reset summation value
    summation = 0;
    for t = 0:Tmax
        summation = summation + abs(I(t+1) - (N*I(1))/(I(1) + (N-I(1))*exp(-1*betahat*t)))^2;
    end

    ideal_I(i,1) = summation;
    i = i + 1;
    N = N + 1;
end

figure(6)
plot(X, ideal_I)
min(ideal_I)

%% Exercise 2 Part 1)
% t0 is all the way from Exercise where t0 = 52 represents the first time
% that the number of infected is greater than or equal to 5.
% vector I will store the 119 days of infection after the day t0
% tau is the incubation/infection period
I = zeros(120,1);
tau = 7;

% Here we are copying over 119 ADDITIONAL days so there's 120 days in total
for j = 0:119
    I(j+1, 1) = x(t0 + j + tau, 1) - x(t0 + j - tau, 1);
end

figure(7)
plot(I)

%% Excercise 2 part 2)
% Here new_N is just N in the context of the problem. I wanted to define a
% new variable because N was used earlier in the homework. 
syms alpha beta new_N



%% Exercise 2 part 3)
a = [1/10 1/9 1/8 1/7 1/6 1/5];
b = [0.8 0.9 0.95 1 1.05 1.1 1.15 1.2 1.3 1.4 1.5 1.6];
c = Nmax * [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];
JABN = zeros(length(a), length(b), length(c));


for i = 1:length(a)
    for j = 1:length(b)
        for m = 1:length(c)
            Isim = SIR_euler(a(i), b(j)*a(i), c(m), I(t0), Tmax);

            summation = 0;
            
            for t = 0:Tmax
                summation = summation + abs(I(t+1) - Isim(t+1))^2;
            end

            JABN(i,j,m) = summation;
        end
    end
end


%%
B = repelem(transpose(b),[1],[10])
C = repelem(c,[12],[1])
surf(B,C, squeeze(JABN(1,:,:)))

B = repelem(transpose(b),[1],[10])
C = repelem(c,[12],[1])
surf(B,C, squeeze(JABN(2,:,:)))
