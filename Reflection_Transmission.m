
% Reflection and Transmission Coefficients at a layer boundary

close all
clear all
im=-sqrt(-1);
alfa1=3.9;
alfa2=4.49;
rho1=2.9;
rho2=3.38;

%alfa1=4.5;
%alfa2=6.48;
beta1=alfa1/sqrt(3);
beta2=alfa2/sqrt(3);
%rho1=2.5;
%rho2=2.81;
as1=1/(alfa1*alfa1);
as2=1/(alfa2*alfa2);
bt1=1/(beta1*beta1);
bt2=1/(beta2*beta2);
pend=1/beta1;
pbeg=0;
np=100;
dp=(pend-pbeg)/(np-1);
mu1=rho1*beta1*beta1;
mu2=rho2*beta2*beta2;
for ip=1:np
  p=pbeg+(ip-1)*dp;
  rp(ip)=p;
  ps=p*p;
  a1=vwnb(as1,ps);
  a2=vwnb(as2,ps);
  b1=vwnb(bt1,ps);
  b2=vwnb(bt2,ps);
  tm1=2*mu1*p;
  tm2=2*mu2*p;
  tm11=tm1*p-rho1;
  tm22=tm2*p-rho2;
  rm=[-p -b1 p -b2; a1 -p a2  p; -tm11 -tm1*b1 tm22 -tm2*b2;...
                  tm1*a1 -tm11 tm2*a2 tm22];
  irm=inv(rm);
%..Incident P wave.....................
  vec1=[p a1 tm11 tm1*a1]';
%..Incident SV wave....................
  vec2=[-b1 p -tm1*b1 tm11]';
  rtip=irm*vec1;
  rtis=irm*vec2;
%...Incident P wave ...................
  rpp(ip)=abs(rtip(1,1));
  rps(ip)=abs(rtip(2,1));
  tpp(ip)=abs(rtip(3,1));
  tps(ip)=abs(rtip(4,1));  
%...Incident P wave ...................
  rsp(ip)=abs(rtis(1,1));
  rss(ip)=abs(rtis(2,1));
  tsp(ip)=abs(rtis(3,1));
  tss(ip)=abs(rtis(4,1));  
  
end
pp=abs(asin(rp*alfa1)*180/pi);
maxz=max([max(rsp),max(rss),max(tps),max(tss)]);
fig=figure;
title('Reflection-Transmision Coefficients for Incident S wave','FontSize',9,'Color','r')
subplot(2,2,1),plot(rp,rsp,'r','LineWidth',2)
title('Rsp','Color','r')
grid on
ylabel('Amplitude','FontSize',8)
axis([0 max(rp) 0 1.1*max(rsp)]);
subplot(2,2,2),plot(rp,rss,'r','LineWidth',2)
title('Rss','Color','r')
grid on
axis([0 max(rp) 0 1.1*max(rss) ]);
subplot(2,2,3),plot(rp,tsp,'r','LineWidth',2)
ylabel('Amplitude','FontSize',8)
axis([0 max(rp) 0 1.1*max(tsp)]);
title('Tsp','Color','r')
grid on
xlabel('Ray parameter','FontSize',8)
subplot(2,2,4),plot(rp,tss,'r','LineWidth',2)
axis([0 max(rp) 0 1.1*max(tss)]);
xlabel('Ray Parameter','FontSize',8)
title('Tss','Color','r')
grid on
fig=figure;
title('Reflection-Transmision Coefficients for Incident P wave','FontSize',9,'Color','r')
subplot(2,2,1),plot(pp,real(rpp),'r','LineWidth',2)
title('Rpp','Color','r')
grid on
ylabel('Amplitude','FontSize',8)
axis([0 max(pp)*0+90 0 1.1*max(rpp)]);
subplot(2,2,2),plot(pp,real(rps),'r','LineWidth',2)
title('Rps','Color','r')
grid on
axis([0 max(pp)*0+90 0 1.1*max(rps)]);
subplot(2,2,3),plot(pp,real(tpp),'r','LineWidth',2)
ylabel('Amplitude','FontSize',8)
axis([0 max(pp)*0+90 0 1.1*max(tpp)]);
title('Tpp','Color','r')
grid on
xlabel('Incidence Angle (deg)','FontSize',8)
subplot(2,2,4),plot(pp,real(tps),'r','LineWidth',2)
axis([0 max(pp)*0+90 0 1.1*max(tps)]);
xlabel('Incidence Angle (deg)','FontSize',8)
title('Tps','Color','r')
grid on
