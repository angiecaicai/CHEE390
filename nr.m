function [ xn ] = nr( fun,xg,xl,xr,tol,Ki,zn )
%% Initializing
i=0;
nrfail=0;
check=1;
h=1e-4;

%% Loop
while tol<check
    i=i+1;
    fp=-zn(1)*(Ki(1)-1)^2/(1+xg*(Ki(1)-1))^2-zn(2)*(Ki(2)-1)^2/(1+xg*(Ki(2)-1))^2;
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

end

