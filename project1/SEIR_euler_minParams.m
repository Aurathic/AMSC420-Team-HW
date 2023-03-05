function [alphaMinInd, betaMinInd, deltaMinInd, NMinInd, paramMinp] = SEIR_euler_minParams(omega2, pSet, J, gammas)
    pLen = length(pSet);
    paramMinp = cell(pLen,5); % Store {alpha, beta, gamma, delta, N}
    
    % For each p value...
    for pInd = 1:pLen
        p = pSet(pInd);
        Jp = J(:,:,:,:,pInd);
        gammasp = gammas(:,:,:,:,pInd);
    
        % Get index of minimum error 
        [M, Ind] = min(Jp, [], "all");
        [alphaMinInd, betaMinInd, deltaMinInd, NMinInd] = ind2sub(size(Jp),Ind);
    
        % Find parameters at that index
        gammaMin = gammasp(Ind);
        omegaTemp = reshape(omega2, [], 4);
        paramMin = num2cell(omegaTemp(Ind,:));
        [alphaMin, betaMin, deltaMin, NMin] = paramMin{:};
        
        % Print result
        fprintf("For p = %d, the parameters which reduce the error are\n" + ...
            "alpha = %.3f, beta = %.3f, gamma = %.3f, delta = %.3f, N = %g\n" + ...
            "with an error of %f\n\n", ...
            p,alphaMin, betaMin, gammaMin, deltaMin, NMin, M);
        % Store minimum parameter values for Exercise 2 Part 3
        paramMinp(pInd,:) = {alphaMin, betaMin, gammaMin, deltaMin, NMin};
    end
end
