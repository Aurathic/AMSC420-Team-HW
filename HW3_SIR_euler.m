function Isim = SIR_euler(alpha, beta, N, initial_I, Tmax)
    % here f1, f2, and f3 are the functions for SIR model
    f1 = 0;
    f2 = 0;
    f3 = 0;
    h = 0.01;
    % The initial values for S(0), I(0), and R(0)
    x1 = N;
    x2 = initial_I;
    x3 = 0;
    Isim = zeros(1,1);
    

    % Here x1 = S, x2 = I, x3 = R
    for t=0:h:Tmax
        f1 = -beta*x1*(x2/N);
        f2 = beta*x1*(x2/N) - alpha*x2;
        f3 = alpha*x2;
        x1 = x1 + f1*h;
        x2 = x2 + f2*h;
        x3 = x3 + f3*h;

        if mod(t,1) == 0
            Isim(t+1,1) = x2;
        end
    end
    
    Isim;
end