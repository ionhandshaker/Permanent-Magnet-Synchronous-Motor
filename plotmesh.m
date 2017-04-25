function plotmesh(feain)
clf
[m,n]=size(feain.object);
ielob=find(feain.el(:,12)==1);
elob=feain.el(ielob,:);
xm=[feain.xn(elob(:,1:3))';feain.xn(elob(:,1))'];
ym=[feain.yn(elob(:,1:3))';feain.yn(elob(:,1))'];
h=line(xm,ym);
fill(xm,ym,'white');
hold on
for i=1:n
    ielob=find(feain.el(:,12)==i+1);
    elob=feain.el(ielob,:);
    xm=[feain.xn(elob(:,1:3))';feain.xn(elob(:,1))'];
    ym=[feain.yn(elob(:,1:3))';feain.yn(elob(:,1))'];
    h=line(xm,ym);
    fill(xm,ym,feain.object(i).color);
    hold on
end
hold off
axis('equal')