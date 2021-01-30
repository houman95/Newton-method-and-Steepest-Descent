function [x_min, f_min, iter] = SD_GSS(f, gf, x0, Stop_tol, GSS_tol,varargin)
    iter = 0;
    x = x0
    while(1)
        iter = iter+1;
        P = -gf(x);
        lnsrch = @(a) f(x + a*P);
        [alpha, N] =  GSS(lnsrch, 0, 10, GSS_tol);
        x_new = x + alpha*P;
        if(norm(x_new - x,2)<= Stop_tol)
            break
        end
        x = x_new;
    end
    x_min = x_new;
    f_min = f(x_min);    
end