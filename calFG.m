function [F,G]=calFG(x,y,mmat)
mu0=4e-7*pi;
N=length(mmat);
F=sparse(N,N);
G=sparse(N,N);
for i=1:N
    m=mmat(i);
    F(i,i)=m/y*((x/y)^m+(y/x)^m)/((x/y)^m-(y/x)^m)/mu0;
    G(i,i)=2*m/y/((y/x)^m-(x/y)^m)/mu0;
end
