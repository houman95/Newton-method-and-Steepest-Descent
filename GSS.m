function [x_min, N] =  GSS(f, a, b, epsilon)
    %% initialization
    l = a;
    h = b;
    %% check the number of required iteration 
    rho = (3 - sqrt(5))/2; %according to the alg
    N = ceil(log(epsilon/(b-a))/log(1-rho)); % # of iterations
    if 1>N % i.e. the initial interval is presice enough
        x_min = 0.5*(l+h);
    else
    %% first step
        t_l = l + (h-l)*rho;
        t_h = h - (h-l)*rho;
        fol = f(t_l);
        foh = f(t_h);
        if fol<foh
            h = t_h;
            
            foh = fol;
            fol = f(l + (h-l)*rho);
            t_h = t_l;
            t_l = l + (h-l)*rho;
        else
            l = t_l;           
            fol = foh;
            foh = f(h-(h-l)*rho);
            t_l = t_h;
            t_h = h-(h-l)*rho;
        end
        if N==1  % only one step is required
            x_min = 0.5*(l+h);
        else   
            %% next steps (if any)
            for i = 2:N
%                 t_l = l + (h-l)*rho;
%                 t_h = h - (h-l)*rho;
                if fol<foh
                    h = t_h;
                    foh = fol;
                    fol = f(l + (h-l)*rho);
                    t_h = t_l;
                    t_l = l + (h-l)*rho;
                else
                    l = t_l;           
                    fol = foh;
                    foh = f(h-(h-l)*rho);
                    t_l = t_h;
                    t_h = h-(h-l)*rho;
                end
            end
            x_min = 0.5*(l+h);
        end
    end
end