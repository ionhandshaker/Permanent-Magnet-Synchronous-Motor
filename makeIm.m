function Im=makeIm(feain)
Mr=5.57e5;
mu0=4e-7*pi;
elements=length(feain.elreal);
Im=zeros(length(feain.xn),1);
Mx=zeros(elements,1);
My=zeros(elements,1);
N=feain.N;
oddpiece=find(feain.elreal(:,12)==2*N+1+2|feain.elreal(:,12)...
    ==2*N+3+2|feain.elreal(:,12)==2*N+5+2);
xcodd=feain.elreal(oddpiece,13);
ycodd=feain.elreal(oddpiece,14);
thetaodd=atan2(ycodd,xcodd);
evenpiece=find(feain.elreal(:,12)==2*N+2+2|feain.elreal(:,12)...
    ==2*N+4+2|feain.elreal(:,12)==2*N+6+2);
xceven=feain.elreal(evenpiece,13);
yceven=feain.elreal(evenpiece,14);
thetaeven=atan2(yceven,xceven);
Mx(oddpiece)=Mr*cos(thetaodd);
My(oddpiece)=Mr*sin(thetaodd);
mm=length(oddpiece);
for i=1:mm
    t=oddpiece(i);
    Im(feain.elreal(t,1))=Im(feain.elreal(t,1))+...
        (My(t)*feain.bv(t,1)-Mx(t)*feain.gv(t,1))/2;
    Im(feain.elreal(t,2))=Im(feain.elreal(t,2))/2+...
        (My(t)*feain.bv(t,2)-Mx(t)*feain.gv(t,2))/2;
    Im(feain.elreal(t,3))=Im(feain.elreal(t,3))+...
        (My(t)*feain.bv(t,3)-Mx(t)*feain.gv(t,3))/2;
end
mm=length(evenpiece);
for i=1:mm
    t=evenpiece(i);
    Im(feain.elreal(t,1))=Im(feain.elreal(t,1))+...
        (My(t)*feain.bv(t,1)-Mx(t)*feain.gv(t,1))/2;
    Im(feain.elreal(t,2))=Im(feain.elreal(t,2))+...
        (My(t)*feain.bv(t,2)-Mx(t)*feain.gv(t,2))/2;
    Im(feain.elreal(t,3))=Im(feain.elreal(t,3))+...
        (My(t)*feain.bv(t,3)-Mx(t)*feain.gv(t,3))/2;
end
mm=length(evenpiece);
for i=1:mm
    t=evenpiece(i);
    Im(feain.elreal(t,1))=Im(feain.elreal(t,1))+...
        (My(t)*feain.bv(t,1)-Mx(t)*feain.gv(t,1))/2;
    Im(feain.elreal(t,2))=Im(feain.elreal(t,2))+...
        (My(t)*feain.bv(t,2)-Mx(t)*feain.gv(t,2))/2;
    Im(feain.elreal(t,3))=Im(feain.elreal(t,3))+...
        (My(t)*feain.bv(t,3)-Mx(t)*feain.gv(t,3))/2;
end

    