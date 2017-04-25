function feaout=torqueplot_real(feain)
tic
feaout=feain;
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
mu0=4e-7*pi;
[D,P,Dw,Pw]=makeDP(feaout,Sthmn);
KinvKsb=Kss\Ksb;
Bs=D*(Kbs*KinvKsb-Ksbb)*P;
Br=D*(Kbr*(Krr\Krb)-Krbb)*P;
[Fs,Gs]=calFG(rr,rs,mmat);
[Fr,Gr]=calFG(rs,rr,mmat);
stf11=Bs-(2*pi*rs/Nt)*Fs;
stf22=Br+(2*pi*rr/Nt)*Fr;
Im=makeIm(feain);
Imr=Im(feain.sb+1:feain.rt);
Imb=Im(feain.rt+1:feain.rb);
rght2=D*Kbr*(Krr\(-Imr))+D*Imb;
M=Ismatrix(feain);
KinvM=Kss\M;
rght1=D*Kbs*KinvM;
clear torque
ri=feaout.rso;
ro=feaout.rri;
i1=1;
nhs=length(mmat);
Nturn=3;
Objslot=[1,4,7,10,13,16;5,2,11,8,17,14;3,6,9,12,15,18];
objstart=19;
for i=1:3
    Sslot=0;
    C1=[];C2=[];
    for j=1:6
        slotelem=find(feain.el(:,12)==objstart+Objslot(i,j));
        mm=length(slotelem);
        Ael=feain.el(slotelem,4);
        if mod(j,2)==0
            xi=-1;
        else
            xi=1;
        end
        C1=[C1,Ael'./3*xi];
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
endangle=500;
for angle=0:5:endangle
    thetax=angle*pi/180;
    T=sparse(nhs,nhs);
    Tinv=sparse(nhs,nhs);
    for i=1:nhs
        m=mmat(i);
        T(i,i)=cos(m*thetax);
        if(mod(i,2)==0)
            T(i,i-1)=-sin(m*thetax);
        else
            T(i,i+1)=sin(m*thetax);
        end
    end
    Tinv=T';
    stf12=-(2*pi*rs/Nt)*Gs*T;
    stf21=(2*pi*rr/Nt)*Gr*Tinv;
    feaout.thetax=thetax;
    Jslot=makeJslot(feaout);
    rght11=rght1*Jslot';
    rght=[rght11;rght2];
    
    Agaptmp=stf\rght;
    
    torq=[];
    Nh=length(mmat);
    As=Agaptmp(1:Nh)/Nt*2;
    As(Nh-1)=As(Nh-1)/2;
    Ascomplex=As(1:2:Nh-1)+1i*As(2:2:Nh);
    Arcomplex=Ar(1:2:Nh-1)+1i*Ar(2:2:Nh);
    Br_ricomplex=-1/ri*1i*mmat(1:2:Nh-1).*Ascomplex;
    Htheta_ri=Fs*As+Gs*T*Ar;
    Btheta_ri=mu0*Htheta_ri;
    Btheta_ricomplex=Btheta_ri(1:2:Nh-1)+1i*Btheta_ri(2:2:Nh);
    Bterm=real(Br_ricomplex.*conj(Btheta_ricomplex));
    torq=-pi*ri^2/mu0*Bterm;
    torque(ii)=sum(torq);
    lambda_a(ii)=lambda1(1,:)*Jslot'-lambda2(1,:)*Agaptmp(1:Nh);
    lambda_b(ii)=lambda1(2,:)*Jslot'-lambda2(2,:)*Agaptmp(1:Nh);
    lambda_c(ii)=lambda1(3,:)*Jslot'-lambda2(3,:)*Agaptmp(1:Nh);
    ii=ii+1;
end
toc
feaout.lambda_a=lambda_a;
feaout.lambda_b=lambda_b;
feaout.lambda_c=lambda_c;
feaout.torque=torque;

    