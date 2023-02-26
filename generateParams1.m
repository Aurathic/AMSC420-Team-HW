function [alphaLen, betaLen, NLen, omega] = generateParams1(Nmax)
    % Compute set omega of possible values for (alpha, beta, N)
    alphaSet = 0.05:0.01:0.4;
    R0Set = 1.5:0.1:1.9;
    NfracSet = 0.02:0.01:0.1;

    alphaLen = length(alphaSet);
    betaLen = length(R0Set);
    NLen = length(NfracSet);

    omega = zeros(alphaLen, betaLen, NLen, 3);
    for alphaInd = 1:alphaLen
        for R0Ind = 1:betaLen
            for NfracInd = 1:NLen
                alpha = alphaSet(alphaInd);
                beta = alpha * R0Set(R0Ind);
                N = Nmax * NfracSet(NfracInd);
                omega(alphaInd, R0Ind, NfracInd, :) = [alpha, beta, N];
            end
        end
    end
end