%% Preprocessing
clear global

% Reading in the Excel file as a table then converting it into a matrix
mytable = readtable('project8_data.xlsx');
x = mytable{2:3,13:end};
x = transpose(x);

% 1st column: cumulative number of detected infections
V = x(:,1)
% 2nd column: Covid-19 related deaths reported in that specific county/city
Y = x(:,2)

%%
Nmax = 236842; % Max population; from Population column
Tmax = 119;    % Number of days we will attempt to model
Vmin = 5;      % See below
tau0 = 7;      % Time between infection and full symptom onset
h = 0.01;      % Step size

% Sets used in omega set generation 
alphaSet = 0.05:0.01:0.2;
R0Set = 1.5:0.1:1.9;
NfracSet = 0.02:0.01:0.1;
% Norm used in error calculation
pSet = [1 2 inf];
pLen = length(pSet);

%%
% Get the first day where at least Vmin were detected as infected
for i = 1:size(x, 1)
    if x(i,1) >= 5
        break
    end
end
t0 = i

% Preprocess rate of infections
I = zeros(Tmax+1,1); % note that I(t) represent the value of I at t+1
for t=0:Tmax
    I(t+1) = V(t+t0+tau0) - V(t+t0-tau0);
end
I;

I0 = I(1)

%%
% Get parameters
% Function defined in generateParams1.m
[alphaLen, betaLen, NLen, omega1] = ...
    generateParams1(Nmax);

%%
paramMinp = cell(pLen,4); % Store {alpha, beta, gamma, N}

% For each p value...
for pInd = 1:pLen
    p = pSet(pInd);
    Jp = J(:,:,:,pInd);
    gammasp = gammas(:,:,:,pInd);

    % Get index of minimum error 
    [M, Ind] = min(Jp, [], "all");
    [alphaMinInd, betaMinInd, NMinInd] = ind2sub(size(Jp),Ind);

    % Find parameters at that index
    gammaMin = gammasp(Ind);
    omegaTemp = reshape(omega1, [], 3);
    paramMin = num2cell(omegaTemp(Ind,:));
    [alphaMin, betaMin, NMin] = paramMin{:};
    
    % Print result
    fprintf("For p = %d, the parameters which reduce the error are\n" + ...
        "alpha = %.3f, beta = %.3f, gamma = %.3f,N = %g\n" + ...
        "with an error of %f\n\n", ...
        p,alphaMin, betaMin, gammaMin, NMin, M);
    % Store minimum parameter values for Exercise 1 Part 3
    paramMinp(pInd,:) = {alphaMin, betaMin, gammaMin, NMin};
end

%% Exercise 1 Part 3
for pInd = 1:pLen
    p = pSet(pInd);
    % Get minimum parameter values
    paramMin = paramMinp(pInd,:);
    [alphaMin, betaMin, gammaMin, NMin] = paramMin{:};
    % Recalculated because it's not stored previously
    [Ssim, Isim, Rsim] = SIR_euler(I0, Tmax, alphaMin, betaMin, NMin);

    % (i) Compare Isim and I
    figure
    plot(Isim)
    hold on
    plot(I)
    hold off
    legend(["predicted", "actual"])
    title(['I_{sim} vs I, p = ', num2str(p)])

    % (ii) Compare Y and Ysim = gamma * Rsim
    figure
    plot(gammaMin * Rsim)
    hold on
    plot(Y(t0:t0+Tmax))
    hold off
    legend(["predicted", "actual"])
    title(['Y_{sim} vs Y, p = ', num2str(p)])
end
