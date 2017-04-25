function y=Solveforvoltage(ik,feaout,thetar)
feain=feaout;
Kss=feaout.Kss;
Ksb=feaout.Ksb;
Kbs=Ksb';
Ksbb=feaout.Ksbb;
Krr=feaout.Krr;
Krb=feaout.Krb;
Kbr=Krb';
Krbb=feaout.Krbb;
rs=feaout.rso;
rr=feaout.rri;

Nt=feaout.sb-feaout.st;
N=Nt;
mmat=zeros(N,1);
mmat(1:2:Nt-1)=1:1:Nt/2;
mmat(2:2:Nt)=1:1:Nt/2;
Sthmn=1:1:N;

[D,P,Dw,Pw]=makeDP(feaout,Sthmn);
KinvKsb=Kss\Ksb;
Bs=D*(Kbs*KinvKsb-Ksbb)*P;
Br=D*(Kbr*(Krr\Krb)-Krbb)*P;
[Fs,Gs]=calFG(rr,rs,mmat);
[Fr,Gr]=calFG(rs,rr,mmat);

stf11=Bs-(2*pi*rs/Nt)*Fs;
stf22=Br+(2*pi*rr/Nt)*Fr;

Im=makeIm(feaout);
Imr=Im(feain.sb+1:feaout.rt);
Imb=Im(feain.rt+1:feain.rb);
rght2=D*Kbr*(Krr\(-Imr))+D*Imb;

M=Ismatrix(feaout);
KinvM=Kss\M;
rght1=D*Kbs*KinvM;

Nturn=3;
Objslot=[1,4,7,10,13,16;5,2,11,8,17,14;3,6,9,12,15,18];
objstart=19;

for i=1:3
    Sslot=0;
    C1=[];C2=[];
    for j=1:6
        slotelem=find(feain.el(:,12)==objstart+Objslot(i,j));
        mm=length(slotelem);
        Ae1=feain.el(slotelem,4);
        if(mod(j,2)==0
            xi=-1;
        else
            xi=1;
        end
        C1=[C1,Ae1'./3*xi];
        C2temp=zeros(mm,length(Kss));
        for t=1:mm
            C2temp(t,feain.el(slotelem(t),1))=1;
            C2temp(t,feain.el(slotelem(t),2))=1;
            C2temp(t,feain.el(slotelem(t),3))=1;
        end
        C2=[C2;C2temp];
        Sslot=Sslot+sum(Ael);
    end
    tt=Nturn/Sslot*C1*C2;
    lambda1(i,:)=tt*KinvM;
    lambda2(i,:)=tt*KinvKsb*P;
end
nhs=length(mmat);
T=sparse(nhs,nhs);
Tinv=sparse(nhs,nhs);
for i=1:nhs
    m=mmat(i);
    T(i,i)=cos(m*thetar);
    if mod(i,2)==0
        T(i,i-1)=-sin(m*thetar);
    else
        T(i,i+1)=sin(m*thetar);
    end
end
Tinv=T';
stf12=-(2*pi*rs/Nt)*Gs*T;
stf21=(2*pi*rr/Nt)*Gr*Tinv;
stf=[stf11 stf12; stf21 stf22];

Jslot=makeJslot4v(feaout,ik);
rght11=rght1*Jslot';
rght=[rght11;rght2];

Agaptmp=stf\rght;
lambda_a=lambda1(1,:)*Jslot'-lambda2(1,:)*Agaptmp(1:nhs);
lambda_b=lambda1(2,:)*Jslot'-lambda2(2,:)*Agaptmp(1:nhs);
lambda_c=lambda1(3,:)*Jslot'-lambda2(3,:)*Agaptmp(1:nhs);

R=9e-7;
f=60;
ts=feaout.ts;
y=[lambda_a;lambda_b;lambda_c];

