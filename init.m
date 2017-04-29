%%% Code for meshing

function feaout=init(feain)
feaout=feain;
disp('Generating initial mesh...');
[m,n]=size(feaout.object);
p=[];
for i=1:n
    p=[p;feaout.object(i).bx feaout.object(i).by];
end
%sort points
[m,n]=size(p);
r=[];
err=1e-12;
for i=1:m
    for j=i+1:m
        if(abs(p(i,:)-p(j,:))<err)
            r=[r j];
        end
    end
end
feaout.xn=p(:,1);
feaout.yn=p(:,2);
feaout.xn(r)=[];
feaout.yn(r)=[];
%%%
xmin=min(p(:,1));
xmax=max(p(:,2));
ymin=min(p(:,2));
ymax=max(p(:,2));
feaout.xn=[xmin-1; xmax+1; xmax+1; xmin-1; feaout.xn];
feaout.yn=[ymin-1; ymin-1; ymax+1; ymax+1; feaout.yn];
feaout.el=zeros(2,14);
feaout.el(1,1:3)=[1 2 3];
feaout.el(2,1:3)=[1 3 4];
feaout.el(1:2,:)=elcalc(feaout.el(1:2,:),...
    feaout.xn(1:4),feaout.yn(1:4),1);
m=length(feaout.xn)-4;
ln=0;
for i=1:m
    ln=ln+1;
    feaout.el=insertnode(feaout.el,feaout.xn,feaout.yn,4+i,1);
end
[i1,j1]=find(feaout.el(:,1:3)==1);
[i2,j2]=find(feaout.el(:,1:3)==2);
[i3,j3]=find(feaout.el(:,1:3)==3);
[i4,j4]=find(feaout.el(:,1:3)==4);
ir=[i1;i2;i3;i4];
feaout.el(ir,:)=[];
feaout.el(:,1:3)=feaout.el(:,1:3)-4;
feaout.xn(1:4)=[];
feaout.yn(1:4)=[];
%%%
feaout.el(:,12)=1;
[me,ne]=size(feaout.el);
[mo,no]=size(feaout.object);
for i=1:me
    for j=1:no
        angl=zeros(1,1);
        N=length(feaout.object(j).bx);
        for k=1:N-1;
            v1=(feaout.object(j).bx(k)-feaout.el(i,13))+1j*(feaout.object(j).by(k)-feaout.el(i,14));
            v2=(feaout.object(j).bx(k+1)-feaout.el(i,13))+1j*(feaout.object(j).by(k+1)-feaout.el(i,14));
            v2p=1/abs(v1)*((real(v1)*real(v2)+imag(v1)*imag(v2)));
            ang2=angle(v2p);
            if(ang2>pi)
                ang2=ang2-2*pi;
            end
            angl=angl+ang2;
        end
        v1=(feaout.object(j).bx(N)-feaout.el(i,13))+1j*(feaout.object(j).by(N)-feaout.el(i,14));
        v2=(feaout.object(j).bx(1)-feaout.el(i,13))+1j*(feaout.object(j).by(1)-feaout.el(i,14));
        v2p=1/abs(v1)*((real(v1)*real(v2)+imag(v1)*imag(v2))+1j*(-imag(v1)*real(v2)+real(v1)*imag(v2)));
        ang2=angle(v2p);
        if(ang2>pi)
            ang2=ang2-2*pi;
        end
        angl=angl+ang2;
        err1=1e-1;
        if(abs(abs(angl)-2*pi)<err1)
            feaout.el(i,12)=j+1;
            break
        end
    end
end


