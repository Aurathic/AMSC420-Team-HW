function objFunc = objectiveFunction(Y, Rsim, I, Isim, gamma, lambda, p)
    %size(Y)
    %size(Rsim)
    %size(I)
    %size(Isim)
    objFunc = norm(I-Isim, p) + lambda * norm(Y - gamma * Rsim, p);
end