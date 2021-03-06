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
z=zeros(1,2);
%prompt=('input the feed fraction of isobutane and carboon dioxide in a 1x2 matrix ');
%z=input(prompt);
z=[0.7 0.3];
zn=z./sum(z);

%% Calculate alpha and ap (a for pure component)
kcst=0.37464+(1.54226.*w)-(0.26992.*(w.^2));
alpha=(1+kcst.*(1-sqrt(Tr))).^2;
ap=ac.*alpha;
bp=bc;
%%
[ xr,yr,vap,liq,fail ] = calflash( 2.923e6,T,Pc,Tc,zn,w,ap,bp,ip )
