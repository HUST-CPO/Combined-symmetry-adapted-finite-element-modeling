function  PlotnormE2(node,elem,Ex,Ey,Ez,alpha)

normE=sqrt(abs(Ex.*Ex)+abs(Ey.*Ey)+abs(Ez.*Ez));
h2=trisurf(elem,node(:,1),node(:,2),normE);
hold on;
trisurf(elem,node(:,1),-node(:,2),normE);
hold on;
shading interp
set(h2,'edgecolor','black','edgealpha',alpha);
view(2); axis equal; axis tight; axis off;

colormap(jet)
% colorbar

caxis([min(normE) max(normE)]);

end