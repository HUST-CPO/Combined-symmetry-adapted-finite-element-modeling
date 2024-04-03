function  PlotE(node,elem,Ez,alpha)
figure(1);
h2=trisurf(elem,node(:,1),node(:,2),Ez);
shading interp
set(h2,'edgecolor','black','edgealpha',alpha);
view(2); axis equal; axis tight; axis off;

colormap jet;
colorbar

minE=min(Ez);
maxE=max(Ez);
if minE==maxE

else
    caxis([minE maxE]);
end

end