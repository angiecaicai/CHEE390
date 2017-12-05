function [zn] = newtonrm(fun, xi, xf, dx0, tol)

xl=xi;
yl=feval(fun, xl);
xr=xl+dx0;
yr=feval(fun, xr);
zn=zeros(1,100);
nrf=0;

while xr<xf
    while sign(yr*yl)==1 && xr<xf
        xl=xr;
        xr=xl+dx0;
        yl=yr;
        yr=feval(fun,xr);
    end
    
    xg=(xl*yr-xr*yl)/(yr-yl);
    [xn,nrfail]=newton(fun,xg,xl,xr,tol);

    if (nrfail==0)
        nrf=nrf+1;
        zn(1,nrf) = xn;
    else
        [xn,sing]=bisection(fun,xl,xr,tol);
        if sing == 0
            nrf = nrf +1;
            zn(1,nrf) = xn;
        end
    end

    xl = xr;
    xr = xl +dx0;
    yl = yr;
    yr = feval(fun,xr);
end
zn(nrf+1:length(zn))=[];

    


