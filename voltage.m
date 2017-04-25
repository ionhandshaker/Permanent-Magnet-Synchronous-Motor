function feaout=voltage(feain)
feaout=feain;
Vmag=110;
R=9e-7;
f=60;
ts=1/f/72;
feaout.ts=ts;
P=3;
thetae0=0;
unit3=[1;1;1];
vk=0*unit3;
ik=0*unit3;
F1k=0*unit3;
thetar0=0;
thetar=thetar0;
y=Solveforvoltage(0*unit3,feaout,thetar);
F2k=y-R*ts/2*[1;1;1];
m=0;
feaout.ik=zeros(3,180);

for t=ts:ts:72*ts
    thetae=thetae0+2*pi*f*t;
    vkp1=Vmag*[sin(thetae);sin(thetae-2*pi/3);sin(thetae+2*pi/3)];
    thetarp1=thetaar0+2*pi/P*f*t;
    F2kp1=Solveforvoltage(0*unit3,feaout,thetarp1)'
    rght=ts/2*(vkp1+vk)-(R*ts/2*[1;1;1]-F1k).*ik+(F2k-F2kp1);
    restart=[];tol=[];maxit=[];M1=[];M2=[];
    x0=[];
    ikp1=gmres(@(x)Solveforvoltage1(x,feaout,thetarp1),rght,restart,tol,maxit);
    y=Solveforvoltage1(ikp1,feaout,thetarp1);
    F1kp1=(y-R*ts*[1;1;1].*ikp1)./ikp1;
    m=m+1;
    feaout.ik(:,m)=ikp1;
    vk=vkp1;
    ik=ikp1;
    F1k=F1kp1;
    F2k=F2kp1;
end

t=5:5:360;
mdl=feaout;
plot(t,mdl.ik(1,1:72),t,mdl.ik(2,1:72),':',t,mdl.ik(3,1:72),'-.');
legend('Phase A','Phase B','Phase C')
xlabel('Electrical Phase in Degree')
ylabel('Amp/m^2')
title('Current Density')
va=110*sin(thetae*pi/180);
vb=110*sin(thetae*pi/180-2*pi/3);
vc=110*sin(thetae*pi/180+2*pi/3);
plot(t,va,t,vb,':',t,vc,'-.')
legend('Phase A','Phase B','Phase C')
xlabel('Electrical Phase in Degree')
ylabel('Volt')
title('Input Voltage')