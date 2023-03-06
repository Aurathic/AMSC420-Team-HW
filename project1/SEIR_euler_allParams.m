function [J,gammas] = SEIR_euler_allParams(t0, Tmax, I0, Y, omega2, pSet)
    % gammas(alphaInd, betaInd, deltaInd, NInd, p) contains 
    % gamma as calculated to minimize p-norm of residual
    % of actual results as compared to euler model 
    % with given parameters (alpha, beta, delta, N)
    % J(...) contains the minimized function value

    % Get number of variations of each parameter
    omegaLens = num2cell(size(omega2));
    [alphaLen, betaLen, deltaLen, NLen, ~] = omegaLens{:};
    pLen = length(pSet);

    gammas = zeros(alphaLen, betaLen, deltaLen, NLen, pLen);
    J = zeros(alphaLen, betaLen, deltaLen, NLen, pLen);
    
    % WARNING: Takes a while to calculate
    fprintf("alphas iterated over: ")
    for alphaInd = 1:alphaLen
        fprintf("%.2f, ", omega2(alphaInd,1,1,1,1))
        for betaInd = 1:betaLen
            for deltaInd = 1:deltaLen
                for NInd = 1:NLen
                    % Get parameters from set, use Euler scheme
                    params = num2cell(squeeze(omega2(alphaInd, betaInd, deltaInd, NInd, :)));
                    [alpha, beta, delta, N] = params{:};
                    [Ssim, Esim, Isim, Rsim] = SEIR_euler(I0, Tmax, alpha, beta, delta, N);
                    
                    % For each p, find gamma minimizing p-norm & store
                    for pInd = 1:pLen
                        p = pSet(pInd);
                        % Function defined in minimizeGamma.m
                        [gamma, ~] = minimizeGamma(t0, Tmax, Y, Rsim, p);
                        gammas(alphaInd, betaInd, deltaInd, NInd, pInd) = gamma;
                        
                        J(alphaInd, betaInd, deltaInd, NInd, pInd) = minVal;
                    end
                end
            end
        end
    end
end

