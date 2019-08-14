function par = kdeo(kp)
%
%   deorientation of kp
%
%deorientation
%step 1: cloude parameters
fee1=angle(kp(1));
fee2=angle(kp(2));
fee3=angle(kp(3));
kpN=abs(kp);
beta=atan(kpN(3)/kpN(2));
alpha=atan( sqrt(kpN(2)^2+kpN(3)^2) / kpN(1) );
%step 2: cross minimization
psim=calpsim(beta,fee2,fee3);
U=[ 1       0           0
    0 cos(2*psim)  sin(2*psim)
    0 -sin(2*psim) cos(2*psim)];
kpd=U*kp;
Upl=[0.5  0.5 0
    0    0   0.707
    0.5 -0.5 0];
kld=Upl*kpd;
%step 3: parameterization
if abs(kld(3))>abs(kld(1))
    z = kld(1)/kld(3);
else
    z = kld(3)/kld(1);
end
a=atan( abs( kld(3)/kld(1) ) );
b=angle( kld(3)/kld(1) )/2;
c=mod(atan( sqrt( abs(kld(1))^2 + abs(kld(3))^2 ) / abs(kld(2)) ), pi);
w=cos(c);
u0=cos(2*a);
v0=cos(2*b);
psi=psim*180/pi;
if u0<0
    psi=psi+90;
end
u0=abs(u0);

if isnan(v0)
    v0=0;
end
if isnan(u0)
    u0=0;
end
par = [real(z),imag(z),psi,u0,v0,w];


function psim=calpsim(beta,fee2,fee3)
%   calculate psim
cdp=cos(fee2-fee3);
tbe=tan(2*beta);
x=2*beta-mod(2*beta,pi)+mod(atan(tan(2*beta)*abs(cdp)),pi);
psim=mod(x*sign(cdp)/4,pi/2);