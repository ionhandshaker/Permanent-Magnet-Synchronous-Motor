%%Based on the work of Danhong, University of Pennsylvania
function feaout=segmdl
xx=0;
mu0=4e-7*pi;
J=[0 -1;1 0];
dtr=pi/180;
unit=0.001;
disp('Generating object...')
N=18;
thetat=360/N;
thetapf=thetat;
rsb=2.99;
rsi=4.23;
rso=6.565;
rg=7.3964;
rm=8.2183;
rbr=10.2;
seg=6;
rsa=rso+3;
%%%%
w1=0.149;
t=0.5;
%%%%
slot(:,1)=[rso;0];
slot(:,2)=[rso-w1;-t/2];
slot(:,3)=[rsi+0.05;-t/2];
slot(:,4)=[rsi;0];
slot(:,5:8)=[fliplr(slot(1,:));fliplr(-slot(2,:))];
%%%
pt0(:,1:4)=slot(:,5:8);
pt0(:,5:24)=rso*[cos(1*dtr:1*dtr:20*dtr);sin(1*dtr:1*dtr:20*dtr)];
pt0(:,25:28)=expm(J*thetat*dtr)*slot(:,1:4);
pt0(:,29:31)=rsi*[cos(thetat*dtr*3/4:-thetat*dtr/4:thetat*dtr/4);
                  sin(thetat*dtr*3/4:-thetat*dtr/4:thetat*dtr/4)];
%%%%
slot=expm(-J*45*dtr)*slot;
pt0=expm(-J*45*dtr)*pt0;

pbs=[];
for i=1:N
    temp1=expm(-J*(i-1)*thetat*dtr)*pt0(:,28:31);
    temp2=expm(-J*(i-1)*thetat*dtr)*pt0(:,1);
    pbs=[pbs temp1 temp2];
end
%%%stator tooth
for i=1:N
    pt=expm(-J*(i-1)*thetat*dtr)*pt0;
    feaout.object(i).name='stator tooth';
    feaout.object(i).bx=pt(1,:)'*unit;
    feaout.object(i).by=pt(2,:)'*unit;
    feaout.object(i).mat=5;
    feaout.object(i).state=1;
    feaout.object(i).color='blue';
end
%%%stator winding slot
for i=1:N
    pslot=expm(-J*(i-1)*thetat*dtr)*slot;
    feaout.object(i+N).name='stator slot';
    feaout.object(i+N).bx=pslot(1,:)'*unit;
    feaout.object(i+N).by=pslot(2,:)'*unit;
    feaout.object(i+N).mat=3;
    feaout.object(i+N).state=1;
    feaout.object(i+N).color='yellow';
end

theta=0:thetat/2:360;
pbo(1,:)=(rsb)*cos(theta*dtr);
pbo(2,:)=(rsb)*sin(theta*dtr);
pb=[pbo pbs];
feaout.object(2*N+1).name='Stator Back Iron';
feaout.object(2*N+1).bx=pb(1,:)'*unit;
feaout.object(2*N+1).by=pb(2,:)'*unit;
feaout.object(2*N+1).mat=5;
feaout.object(2*N+1).state=1;
feaout.object(2*N+1).color='blue';
%%%
startangle1=0;
startangle2=59;
for i=1:seg
    thetaangle(i,:)=(startangle1+(i-1)*60):1:(startangle2+(i-1)*60);
end
for i=1:seg
    segpiece(i).pmi(1,:)=rg*cos(thetaangle(i,:)*dtr);
    segpiece(i).pmi(2,:)=rg*sin(thetaangle(i,:)*dtr);
    segpiece(i).pmo(1,:)=rm*cos(fliplr(thetaangle(i,:)*dtr));
    segpiece(i).pmo(2,:)=rm*sin(fliplr(thetaangle(i,:)*dtr));
    segpiece(i).pm=[segpiece(i).pmo segpiece(i).pmi];
    feaout.object(2*N+1+i).name='Magnet Ring';
    feaout.object(2*N+1+i).bx=segpiece(i).pm(1,:)'*unit;
    feaout.object(2*N+1+i).by=segpiece(i).pm(2,:)'*unit;
    feaout.object(2*N+1+i).mat=4;
    feaout.object(2*N+1+i).state=2;
    feaout.object(2*N+1+i).color='m';
end
theta=0:thetat/2:360;
pbri(1,:)=rm*cos(-theta*dtr)+xx;
pbri(2,:)=rm*sin(-theta*dtr);
pbro(1,:)=rbr*cos(theta*dtr)+xx;
pbro(2,:)=rbr*sin(theta*dtr);
pbr=[pbro pbri];
feaout.object(2*N+seg+2).name='Rotor Back-Iron';
feaout.object(2*N+seg+2).bx=pbr(1,:)'*unit;
feaout.object(2*N+seg+2).by=pbr(2,:)'*unit;
feaout.object(2*N+seg+2).mat=2;
feaout.object(2*N+seg+2).state=2;
feaout.object(2*N+seg+2).color='blue';

pair(1,:)=(rsb)*cos(theta*dtr);
pair(2,:)=(rsb)*sin(theta*dtr);
feaout.object(2*N+seg+3).name='air hole';
feaout.object(2*N+seg+3).bx=pair(1,:)'*unit;
feaout.object(2*N+seg+3).by=pair(2,:)'*unit;
feaout.object(2*N+seg+3).mat=1;
feaout.object(2*N+seg+3).state=1;
feaout.object(2*N+seg+3).color='white';

startangle1=59;
startangle2=60;
for i=1:seg
    gapangle(i,:)=(startangle1+(i-1)*60):1:(startangle2+(i-1)*60);
end
for i=1:seg
    seggap(i).pmai(1,:)=rg*cos(gapangle(i,:)*dtr);
    seggap(i).pmai(2,:)=rg*sin(gapangle(i,:)*dtr);
    seggap(i).pmao(1,:)=rm*cos(fliplr(gapangle(i,:)*dtr));
    seggap(i).pmao(2,:)=rm*sin(fliplr(gapangle(i,:)*dtr));
    seggap(i).pma=[seggap(i).pmai seggap(i).pmao];
    feaout.object(2*N+seg+3+i).name='Magnet Airgap';
    feaout.object(2*N+seg+3+i).bx=seggap(i).pma(1,:)'*unit;
    feaout.object(2*N+seg+3+i).by=seggap(i).pma(2,:)'*unit;
    feaout.object(2*N+seg+3+i).mat=1;
    feaout.object(2*N+seg+3+i).state=2;
    feaout.object(2*N+seg+3+i).color='white';
end
%%%
pap=[];
for i=1:seg
    pap=[pap segpiece(i).pmi seggap(i).pmai];
end
theta=0:1:360;
psair(1,:)=(rsa)*cos(theta*dtr);
psair(2,:)=(rsa)*sin(theta*dtr);
pas=[psair pap];
feaout.object(2*N+2*seg+4).name='air gap';
feaout.object(2*N+2*seg+4).bx=pas(1,:)'*unit;
feaout.object(2*N+2*seg+4).by=pas(2,:)'*unit;
feaout.object(2*N+2*seg+4).mat=1;
feaout.object(2*N+2*seg+4).state=1;
feaout.object(2*N+2*seg+4).color='white';
pa1(:,1:22)=pt0(:,4:25);
for i=1:N
    pa=expm(J*(i-1)*thetat*dtr)*pa1;
    paline(:,22*(i-1)+1:22*i)=pa(:,1:22);
end
pass=[paline psair];
feaout.object(2*N+2*seg+5).name='air gap';
feaout.object(2*N+2*seg+5).bx=pass(1,:)'*unit;
feaout.object(2*N+2*seg+5).by=pass(2,:)'*unit;
feaout.object(2*N+2*seg+5).mat=1;
feaout.object(2*N+2*seg+5).state=1;
feaout.object(2*N+2*seg+5).color='white';
%%%
feaout.mat.number(1)=3;
feaout.mat.mat(1)=2;
feaout.mat.mat(2)=4;
feaout.mat.mat(3)=5;
feaout.mat.mag(1,:)='mofb1';
feaout.mat.mag(2,:)='mofb2';
feaout.mat.mag(3,:)='mofb3';
feaout.mat.dmag(1,:)='dmdb1';
feaout.mat.dmag(2,:)='dmdb2';
feaout.mat.dmag(3,:)='dmdb3';
feaout.mat.dmb(1,:)='dmb1';
feaout.mat.dmb(2,:)='dmb2';
feaout.mat.dmb(3,:)='dmb1';

feaout.rro=rbr*unit;
feaout.rri=rg*unit;
feaout.rso=rsa*unit;
feaout.N=N;
plotobjects(feaout)



