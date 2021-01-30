function [y varargout] = Powel(x)
nout = max(nargout,1)-1;

    persistent PowelCount
    if isempty(PowelCount)
        PowelCount = 0;
    end
    y = (x(1)+10*x(2))^2+5*(x(3)-x(4))^2+(x(2)-2*x(3))^4+10*(x(1)-x(4))^4;
    PowelCount = PowelCount+1;
    global NPowelEval;
    NPowelEval = PowelCount;
    for k=1:nout
        varargout(k) = {PowelCount};
    end
end