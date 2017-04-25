function [D,P,Dw,Pw]=makeDP(feain,Sthmn)
feaout=feain;
Nts=feaout.sb-feaout.st;
for i=1:Nts/2
    for j=1:Nts
        Dw(2*i-1,j)=cos(2*pi*i*(j-1)/Nts);
        Dw(2*i,j)=-sin(2*pi*i*(j-1)/Nts);
    end
end
Pw=2/Nts*Dw';
Pw(:,Nts-1)=1/2*Pw(:,Nts-1);
D=Dw(Sthmn,:);
P=Pw(:,Sthmn);
