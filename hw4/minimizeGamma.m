function [gamma, minVal] = minimizeGamma(t0, Tmax, Y, Rsim, p)
    %MINIMIZEGAMMA Find gamma which minimizes residuals
    % Find argmin for gamma over the p-norm of (Y - gamma*Rsim)
    % p should be 1, 2, or infinity (inf)
    % Also returns the minimized value.
    
    Ycut = Y(t0:t0+Tmax);
    func = @(g) norm(Ycut-g*Rsim, p);
    [gamma, minVal] = fminsearch(func, 1);
end