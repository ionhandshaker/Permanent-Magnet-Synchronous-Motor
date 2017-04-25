function elout=insertnode(el,xn,yn,i,mat)
err=1e-12;
[i1,j1]=find(((xn(i)-el(:,9)).^2+(yn(i)-el(:,10)).^2)<...
    el(:,11)+err);
el1=el(i1,:);
[i2,j2]=find(el1(:,12)==mat);
n=i1(i2);

edges=[el(n,1) el(n,2); el(n,2) el(n,3); el(n,3) el(n,1)];
edges=sort(edges')';
el(n,:)=[];
[m,n]=size(edges);
erm=[];
for j=1:m-1
    for k=j+1:m
        if(edges(j,:)==edges(k,:))
            erm=[erm j k];
        end
    end
end
edges(erm,:)=[];
[med,ned]=size(edges);
newel=zeros(med,14);
newel(:,1:3)=[edges i*ones(med,1)];
newel=elcalc(newel,xn,yn,mat);
elout=[el;newel];