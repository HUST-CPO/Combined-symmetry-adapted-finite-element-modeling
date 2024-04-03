function PlotTri(node,elem)

figure(3);
subplot(1,1,1);
h=trisurf(elem,node(:,1),node(:,2),zeros(size(node,1),1));
set(h,'facecolor',[1 1 1],'edgecolor','black');
view(2); axis equal; axis tight; axis off;

end