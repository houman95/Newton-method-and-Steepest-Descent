function [y varargout] = GR_Rosenbrock(x)
nout = max(nargout,1)-1;

    persistent RosenCount
    if isempty(RosenCount)
        RosenCount = 0;
    end
    y = [-400*(x(1)*x(2) - x(1)^3) + 2*(-1 + x(1));200*(x(2) - x(1)^2)];
    RosenCount = RosenCount+1;
    global NGR_RosenEval;
    NGR_RosenEval = RosenCount;
    for k=1:nout
        varargout(k) = {RosenCount};
    end
end