function [alphaLen, betaLen, deltaLen, NLen, omega] = generateParams2(alphaSet, deltaSet, R0Set, NfracSet, Nmax)
    % Compute set omega of possible values for (alpha, beta, delta, N)
    alphaLen = length(alphaSet);
    betaLen = length(R0Set);
    deltaLen = length(deltaSet);
    NLen = length(NfracSet);

    omega = zeros(alphaLen, betaLen, deltaLen, NLen, 4);
    for alphaInd = 1:alphaLen
        for R0Ind = 1:betaLen
            for deltaInd = 1:deltaLen
                for NfracInd = 1:NLen
                    alpha = alphaSet(alphaInd);
                    beta = alpha * R0Set(R0Ind);
                    delta = deltaSet(deltaInd);
                    N = Nmax * NfracSet(NfracInd);
                    omega(alphaInd, R0Ind, deltaInd, NfracInd, :) = [alpha, beta, delta, N];
                end
            end
        end
    end
end