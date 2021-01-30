function [xmin, fmin, iter] = Newton_GSS(f,gf,Hf,x0,Stop_tol,GSS_tol,varargin)
    iter = 0;
    x = x0;
    while(1)
        iter = iter + 1;
        H = Hf(x);
        try(chol(H));
            
           % fprintf("tried %d",iter);
        catch
            mu=abs(min(eig(H)));
            H = H + (1.0001)*mu*eye(size(H));
            %fprintf("catched %d",iter);
        end        
        p = -inv(H)*gf(x);
        phi = @(a) f(x + a*p);
		derv_phi = @(a)(P'*gf(x+a*P));
		alpha =  bracketingphase(phi,derv_phi, 0.01, 0.9,1,8,GSS_tol);
        x_new = x + alpha*p;
        if(norm(x_new - x,2)<=Stop_tol)
            break
        end
        if(iter > 100)
            break
        end
        x= x_new;
    end
    xmin = x_new;
    fmin = f(xmin);
end