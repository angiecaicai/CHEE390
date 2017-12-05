clear
clc

%% This part can be modified by user base on different components. Here, 1 denotes isobutane while 2 carbon dioxide.
% The unit of temperature is K and the unit of pressure is in Pa.
%The thermodynamics property data is from Introduction to Chemical
%Engineering Thermodynamics.

Tc=[408.1,304.2];
Pc=[36.48*1e5,73.83*1e5];
w=[0.181,0.224];
Zc=[0.282,0.274];

%ip is the interaction parameters given
ip=[0 0.130;0.130 0];

%%Calculate ac, bc (@critical temperature)
ac=(0.45724*(8.314^2)).*(Tc.^2)./Pc;
bc=((0.07780*8.314).*Tc)./Pc;

%% Input T
T=310.928;
Tr=T./Tc;

%% Input feed
zn=[0.7 0.3];

%% Calculate alpha and ap (a for pure component)
kcst=0.37464+(1.54226.*w)-(0.26992.*(w.^2));
alpha=(1+kcst.*(1-sqrt(Tr))).^2;
ap=ac.*alpha;
bp=bc;

%%
deltaP=6894.76;
xresult=zeros(1,100000);
yresult=zeros(1,100000);
Presult=zeros(1,100000);

P=689476;
count=1;
while zn(2)<=0.9
    [ xn,yn,vap,liq,fail ] = calflash( P,T,Pc,Tc,zn,w,ap,bp,ip );
    if vap~=1
        if liq==1 || fail==1
            zn(2)=yn(2);
            zn(1)=1-zn(2);
        elseif vap==0 && liq==0
            xresult(count)=xn(2);
            yresult(count)=yn(2);
            Presult(count)=P/6894.76;
            count=count+1;
        end
    end
    
    P=P+deltaP;
    
end   
    xresult=[0,xresult];
    yresult=[0,yresult];
    Presult=[73.09902,Presult];
    
    xresult(count+1:end)=[];
    yresult(count+1:end)=[];
    Presult(count+1:end)=[];
    plot(xresult, Presult, yresult, Presult);
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
