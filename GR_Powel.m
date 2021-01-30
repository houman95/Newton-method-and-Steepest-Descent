function [y varargout] = GR_Powel(x)
nout = max(nargout,1)-1;

    persistent PowelCount
    if isempty(PowelCount)
        PowelCount = 0;
    end
    y = [2*(x(1)+10*x(2))+40*(x(1)-x(4))^3;20*(x(1)+10*x(2))+4*(x(2)-2*x(3))^3;10*(x(3)-x(4))-8*(x(2)-2*x(3))^3;-10*(x(3)-x(4))-40*(x(1)-x(4))^3];
    PowelCount = PowelCount+1;
    global NGR_PowelEval;
    NGR_PowelEval = PowelCount;
    for k=1:nout
        varargout(k) = {PowelCount};
    end
end
