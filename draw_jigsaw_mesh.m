function draw_jigsaw_mesh(geom,flips,fc,ec,view_angle)

geom.point.coord = geom.point.coord(:,flips) ;
geom.point.coord = geom.point.coord-...
    repmat(mean(geom.point.coord,1),size(geom.point.coord,1),1);  % center the mesh location

pp = geom.point.coord(:,1:3);
t3 = geom.tria3.index(:,1:3);

figure;%('position',[683,-104,1031,676]);
patch('faces',t3,'vertices',pp, ...
    'facecolor',[1 1 1],...
    'edgecolor',ec,...
    'facealpha',1,...
    'linewidth',0.5);

axis equal
axis off
set(gcf,'color','w');
view(view_angle)
