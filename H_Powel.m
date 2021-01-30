function [y varargout] = H_Powel(x)
nout = max(nargout,1)-1;

    persistent PowelCount
    if isempty(PowelCount)
        PowelCount = 0;
    end
    y = [2+120*(x(1)-x(4))^2, 20, 0, -120*(x(1)-x(4))^2
         20, 200+12*(x(2)-2*x(3))^2, -24*(x(2)-2*x(3))^2, 0
         0, -24*(x(2)-2*x(3))^2, 10+48*(x(2)-2*x(3))^2, -10
         -120*(x(1)-x(4))^2, 0, -10, 10+120*(x(1)-x(4))^2];

    PowelCount = PowelCount+1;
    global N_H_PowelEval;
    N_H_PowelEval = PowelCount;
    for k=1:nout
        varargout(k) = {PowelCount};
    end
end
