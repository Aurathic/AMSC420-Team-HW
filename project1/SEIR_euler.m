function [Sout, Eout, Iout, Rout] = SEIR_euler(I0, Tmax, alpha, beta, delta, N)
    h = 0.01; % step size; I assume it's measured in terms of (fractions of) days
    Nsteps = Tmax/h+1; % number of iterations

    % vectors storing the values at each step
    % the ith entry in each vector is the value at the (i-1)th time step,
    % e.g. Ssim(i) is the approximate value of S at time t = h*(i-1)
    Ssim = zeros(Nsteps,1);
    Esim = zeros(Nsteps,1);
    Isim = zeros(Nsteps,1);
    Rsim = zeros(Nsteps,1);
    
    % Initialization
    Ssim(1) = N;
    Esim(1) = I0;
    Isim(1) = I0;
    Rsim(1) = 0;

    Ntot = Ssim(1)+Esim(1)+Isim(1)+Rsim(1);

    % Steps
    for i = 2:Nsteps
        % t = h*(i-1);
        [Sprev, Eprev, Iprev, Rprev] = deal(Ssim(i-1), Esim(i-1), Isim(i-1), Rsim(i-1));
        x = beta * Sprev * Iprev / Ntot;

        dSdt = -x;
        dEdt =  x - delta * Eprev;
        dIdt =  delta * Eprev - alpha * Iprev;
        dRdt =  alpha * Iprev;
        %[dSdt, dEdt, dIdt, dRdt]

        Ssim(i) = Sprev + dSdt * h;
        Esim(i) = Eprev + dEdt * h;
        Isim(i) = Iprev + dIdt * h;
        Rsim(i) = Rprev + dRdt * h;
    end
    % We only care about the values at the day intervals 
    % (every 1/h time-steps)
    Sout = Ssim(1:(1/h):end);
    Eout = Esim(1:(1/h):end);
    Iout = Isim(1:(1/h):end);
    Rout = Rsim(1:(1/h):end);

end