function [y varargout] = Rosenbrock(x)
    nout = max(nargout,1)-1;
    persistent RosenCount
    if isempty(RosenCount)
        RosenCount = 0;
    end
    y = 100*(x(2) - x(1)^2)^2 + (1 - x(1)^2)^2;
    RosenCount = RosenCount+1;
    global NRosenEval;
    for k=1:nout
        varargout(k) = {RosenCount};
    end
end