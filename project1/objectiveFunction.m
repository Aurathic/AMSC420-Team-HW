function objFunc = objectiveFunction(Y, Rsim, I, Isim, gamma, lambda, p)
    %size(Y)
    %size(Rsim)
    %size(I)
    %size(Isim)
    objFunc = norm(I-Isim, p) + lambda * norm(Y - gamma * Rsim, p);
    %if p==inf
    %   objFunc = max(abs(I-Isim)) + ...
    %   lambda * max(abs(Y-gamma*Rsim));
    %else
    %    objFunc = sum(abs(I-Isim).^p) + ... 
    %    lambda * sum(abs(Y-gamma*Rsim).^p);
    %end
end