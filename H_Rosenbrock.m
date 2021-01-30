function [y varargout] = H_Rosenbrock(x)
nout = max(nargout,1)-1;

    persistent RosenCount
    if isempty(RosenCount)
        RosenCount = 0;
    end
    y = [1200*x(1)^2-400*x(2)+2, -400*x(1)
         -400*x(1), 200];
    RosenCount = RosenCount+1;
    global N_H_RosenEval;
    N_H_RosenEval = RosenCount;
    for k=1:nout
        varargout(k) = {RosenCount};
    end
end