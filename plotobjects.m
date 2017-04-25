function plotobjects(feaout)
figure(1)
clf
hold on
[m,n]=size(feaout.object);
for i=1:n
    xm=feaout.object(i).bx/0.001;
    ym=feaout.object(i).by/0.001;
    fill(xm,ym,feaout.object(i).color)
end
axis('equal')