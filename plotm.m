function plotm(v,feain)
clf
[m,n]=size(v);
if(m~=length(feain.el))
    v=v';
end
xm=zeros(3,length(feain.el(:,1)));
ym=zeros(3,length(feain.el(:,1)));
disp('Plotting...')
xm(1,:)=feain.xn(feain.el(:,1))';
ym(1,:)=feain.yn(feain.el(:,1))';
xm(2,:)=feain.xn(feain.el(:,2))';
ym(2,:)=feain.yn(feain.el(:,2))';
xm(3,:)=feain.xn(feain.el(:,3))';
ym(3,:)=feain.yn(feain.el(:,3))';
fill(xm,ym,v')
colorbar
axis('equal')
grid
xlabel('x ','Fontsize',14)
ylabel('y ','Fontsize',14)