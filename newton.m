
function [xn, nrfail] = newton(fun,xg,xl,xr,tol)
% Newton Raphson method for calculating root in the
% interval xrg=[xl xr]
% (c) 2012 Phillip Servio
%% Initializing
i=0;
nrfail=0;
check=1;
h=1e-4;

%% Loop
while tol<check
    i=i+1;
    fp=der(fun,xg,h);
    y=feval(fun,xg);
    xn=xg-y/fp;
    if xn<xl || xn>xr || i>10 || fp==0
        nrfail=1;
        break 
    end
    if abs(xn) <= 1
        check=abs(xg-xn);
    else
        check=abs(1-xg/xn);
    end
    xg=xn;
end