function feaout=makemat(feain)
clear nu
feaout=feain;
mu0=4e-7*pi;
nu0=1/mu0;
nu(1)=nu0;
nu(2)=nu0/20000;
nu(3)=nu0;
nu(4)=nu0;
nu(5)=nu0/20000;
[m,n]=size(feaout.object);
nnodes=length(feaout.xn);
elements=length(feaout.el);
st=feaout.st;
sb=feaout.sb;
rt=feaout.rt;
rb=feaout.rb;
nus=feaout.nus;
biggama=sparse(nnodes,elements);
bigbeta=sparse(nnodes,elements);
clear K
K=sparse(nnodes,nnodes);
for i=1:elements
    if(feaout.el(i,12)~=1)
        feaout.el(i,16)=feaout.object(feaout.el(i,12)-1).mat;
    else
        feaout.el(i,16)=1;
    end
end
mgapc=sqrt(feaout.el(:,13).^2+feaout.el(:,14).^2);
nogap=find(mgapc<feaout.rso|mgapc>feaout.rri);
elreal=feaouut.el(nogap,:);
betav(:,1)=feaout.yn(elreal(:,2))-feaout.yn(elreal(:,3));
betav(:,2)=feaout.yn(elreal(:,3))-feaout.yn(elreal(:,1));
betav(:,3)=feaout.yn(elreal(:,1))-feaout.yn(elreal(:,2));
gammav(:,1)=feaout.xn(elreal(:,3))-feaout.xn(elreal(:,2));
gammav(:,2)=feaout.xn(elreal(:,1))-feaout.xn(elreal(:,3));
gammav(:,3)=feaout.xn(elreal(:,2))-feaout.xn(elreal(:,1));
Ae=elreal(:,4);

clear i j
for i=1:3
    for j=1:3
        K=K+sparse(elreal(:,i),elreal(:,j),nu(elreal(:,16)).*...
            ((betav(:,i).*betav(:,j)+gammav(:,i).*gammav(:,j))./Ae/4)',...
            nnodes,nnodes);
    end
end
Kss=K(1:st,1:st);
Ksb=K(1:st,st+1:sb);
Ksbb=K(st+1:sb,st+1:sb);
Krr=K(sb+1:rt,sb+1:rt);
Krb=K(sb+1:rt,rt+1:rb);
Krbb=K(rt+1:rb,rt+1:rb);

feaout.Kss=Kss;
feaout.Ksb=Ksb;
feaout.Ksbb=Ksbb;
feaout.Krr=Krr;
feaout.Krb=Krb;
feaout.Krbb=Krbb;
feaout.bv=betav;
feaout.gv=gammav;
feaout.K=K(1:nnodes-feaout.bounodenum,1:nnodes-feaout.bounodenum);
feaout.elreal=elreal;
