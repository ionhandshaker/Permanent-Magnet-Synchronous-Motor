function M=Ismatrix(feain)
N=18;
objstart=19;
Jpk=1e7;
st=feaout.st;
M=sparse(st,N);
for j=1:N
    slotelem=find(feain.el(:,12)==objstart+j);
    Ael=feain.el(slotelem,4);
    mm=length(slotelem);
    for i=1:mm
        M(feain.el(slotelem(i),1),j)=M(feain.el(slotelem(i),1),j)+Ael(i);
        M(feain.el(slotelem(i),2),j)=M(feain.el(slotelem(i),2),j)+Ael(i);
        M(feain.el(slotelem(i),3),j)=M(feain.el(slotelem(i),3),j)+Ael(i);
    end
end
M=1/3*M;
