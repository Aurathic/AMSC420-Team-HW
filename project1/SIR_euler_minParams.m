function [alphaMinInd, betaMinInd, NMinInd, paramMinp] = SIR_euler_minParams(omega1, pSet, J, gammas)
    pLen = length(pSet);
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
        % Store minimum parameter values for next part of exercise
        paramMinp(pInd,:) = {alphaMin, betaMin, gammaMin, NMin};
    end
end

