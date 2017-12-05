function [x, y, success] = flash(temperature, pressure, fraction)

N = 2;
R = 8.314;
success = 0;
z=zeros(1,2);
x=zeros(1,2);
y=zeros(1,2);
% initialize Tc Pc omega of 2 species
% 1 is hydrogen and 2 is propane
Tc = zeros(1,2);
Pc = zeros(1,2);
omega = zeros(1,2);
Tc(1) = 33.20;
Tc(2) = 370;
Pc(1) = 1.3e6;
Pc(2) = 4.21e6;
omega(1) = -0.215;
omega(2) = 0.153;
z(1) = fraction;
z(2) = 1-fraction;

m = zeros(1,2);
b = zeros(1,2);
alpha = zeros(1,2);
a = zeros(1,2);
Tr = zeros(1,2);

for i = 1:N
    m(i) = 0.480 + 1.574*omega(i) - 0.176*omega(i)^2;
    Tr(i) = temperature/Tc(i);
    alpha(i) = (1 + m(i)*(1-Tr(i)^0.5))^2;
    a(i) = 0.42747*alpha(i)*R^2*Tc(i)^2/Pc(i);
    b(i) = 0.08664*R*Tc(i)/Pc(i);
end


K = zeros(1,2);
for i = 1:N
    K(i) = exp(log(Pc(i)/pressure)+5.37*(1+omega(i))*(1-Tc(i)/temperature));
end
k=1;

theta = zeros(1,100);
while k<100
    Pzero = z(1)*(K(1)-1)+z(2)*(K(2)-1);
    Pone = 1 - z(1)/K(1) - z(2)/K(2);

    if Pzero*Pone <= 0
        fun=@(g) P_alpha_fun(z,K,g,N);
        solutions = newtonrm(fun, 0, 1, 0.001, 1e-12);
        Nphase = 2;
        alphaP = solutions(1);
    else
        if Pzero <= 0
            alphaP = 0;
            Nphase = 1;
        else
            if Pone >= 0
                alphaP = 1;
                Nphase = 1;
            else
                disp('Error Occurred');
                break;
            end
        end
    end
    
    for i = 1:N
        x(i) = z(i)/(1+alphaP*(K(i)-1));
        y(i) = K(i)*x(i);
    end

    x_sum = sum(x);
    y_sum = sum(y);
    for i = 1:N
        x(i) = x(i)/x_sum;
        y(i) = y(i)/y_sum;
    end

    % av, al, bv and bl
    bv = 0;
    bl = 0;
    av = 0;
    al = 0;
    for i = 1:N
        bv = bv+b(i)*y(i);
        bl = bl+b(i)*x(i);
        for j = 1:N
            av = av + y(i)*y(j)*a(i)^0.5*a(j)^0.5;
            al = al + x(i)*x(j)*a(i)^0.5*a(j)^0.5;
        end
    end

    % calculate A and B and Z for vapor and liquid phases
    Al = pressure*al/(R^2*temperature^2); 
    Bl = pressure*bl/(R*temperature);
    Av = pressure*av/(R^2*temperature^2); 
    Bv = pressure*bv/(R*temperature);
    
    fun_l = @(Z_1) polynomial(Z_1,Al,Bl);
    rl = newtonrm(fun_l, 0, 2, 0.05, 1e-8);
    Zl = min(rl);
    Vl = Zl*R*temperature/pressure;
    
    fun_v = @(Z_2) polynomial(Z_2, Av, Bv); 
    rv = newtonrm(fun_v, 0, 2, 0.05, 1e-8);
    Zv = max(rv);
    Vv = Zv*R*temperature/pressure;

    % fugacity of 2 species at 2 phase
    fl = zeros(1,2);
    fv = zeros(1,2);
    for i=1:N
        fl(i) = pressure*x(i)*exp(b(i)/bl*(Zl-1)-log(Zl-Bl)-Al/Bl*(2*a(i)^0.5/al^0.5-b(i)/bl)*log(1+Bl/Zl));
        fv(i) = pressure*y(i)*exp(b(i)/bv*(Zv-1)-log(Zv-Bv)-Av/Bv*(2*a(i)^0.5/av^0.5-b(i)/bv)*log(1+Bv/Zv));
    end

    for i = 1:N
        theta(k) = theta(k) + y(i)*log(fv(i)/fl(i));
    end
    
    if abs(theta(k)) <= 1e-7
        fprintf('Flash Converged\n');
        success = 1;
        break;
    else
         if k == 1
             for i = 1:N
                 K(i) = K(i)*fl(i)/fv(i);
             end
             if Nphase == 2
                 if k>=100
                     break;
                 else
                     k=k+1;
                 end
             else
                 for i = 1:N
                     K(i)=K(i)*exp(theta(k));
                 end
                 k = k+1;
             end
         else 
             if abs(theta(k)-theta(k-1))<=1e-10
                 if theta(k)>0
                     if alphaP==0
                        fprintf('Feed is Liquid Like\n')
                        break; 
                     end
                 else
                     if alphaP==1
                        fprintf('Feed is Vapor Like\n')
                        break;
                     end
                 end        
             else 
                 for i = 1:N
                     K(i) = K(i)*fl(i)/fv(i);
                 end
                 if Nphase == 2
                         k=k+1;
                 else
                     for i = 1:N
                         K(i)=K(i)*exp(theta(k));
                     end
                     k= k+1;
                 end        
             end
         end
     end
 end



end


