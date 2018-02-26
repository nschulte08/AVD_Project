function [x1,S_x1] = GoldenSection(LB,UB,fun)

LB = 3;
UB = 4;

% alpha = .4;
% PRc = 16;
% Alt = 35000;
% M0 = 1.6;
% AB = 1;

%fun = @(PRf) OptiS(alpha,PRf,PRc,Alt,M0,AB);
x2 = (UB-LB)*.618+LB;
x1 = (UB-LB)*(1-.618)+LB;

while abs(x1-x2)>.01
    x2 = (UB-LB)*.618+LB;
    x1 = (UB-LB)*(1-.618)+LB;
    S_x2 = fun(x2);
    S_x1 = fun(x1);

    if S_x1>S_x2
        LB = x1;
    else
        UB = x2;
    end
    
end

end
