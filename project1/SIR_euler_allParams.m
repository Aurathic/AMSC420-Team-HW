function [J,gammas] = SIR_euler_allParams(t0, Tmax, I0, Y, I, omega1, pSet)
    % gammas(alphaInd, betaInd, NInd, p) contains 
    % gamma as calculated to minimize p-norm of residual
    % of actual results as compared to euler model 
    % with given parameters (alpha, beta, N)
    % J(...) contains the minimized function value
    % WARNING: Takes a while to calculate
    
    % Get number of variations of each parameter
    omegaLens = num2cell(size(omega1));
    [alphaLen, betaLen, NLen, ~] = omegaLens{:};
    pLen = length(pSet);

    gammas = zeros(alphaLen, betaLen, NLen, pLen);
    J = zeros(alphaLen, betaLen, NLen, pLen);
    
    fprintf("alphas iterated over: ")
    for alphaInd = 1:alphaLen
        fprintf("%.2f, ", omega1(alphaInd,1,1,1))
        for betaInd = 1:betaLen
            for NInd = 1:NLen
                % Get parameters from set, use Euler scheme
                params = num2cell(squeeze(omega1(alphaInd, betaInd, NInd, :)));
                [alpha, beta, N] = params{:};
                [Ssim, Isim, Rsim] = SIR_euler(I0, Tmax, alpha, beta, N);
                
                % For each p, find gamma minimizing p-norm & store
                for pInd = 1:pLen
                    p = pSet(pInd);
                    % Function defined in minimizeGamma.m
                    [gamma, ~] = minimizeGamma(t0, Tmax, Y, Rsim, p);
                    gammas(alphaInd, betaInd, NInd, pInd) = gamma;
                    J(alphaInd, betaInd, NInd, pInd) = objectiveFunction(Y(t0:t0+Tmax), Rsim, I(1:Tmax+1), Isim, gamma, 1, p);
                end
            end
        end
    end
end

