function plot_model_shape_src(h,outline,src_loc,tongue_loc,axlim,vh)

% convert to [cm]
outline = outline*1e2;
src_loc = src_loc*1e2;
tongue_loc = tongue_loc*1e2;
axlim = axlim*1e2;

% plot
axes(h)
plot3(outline(:,1),outline(:,2),outline(:,3),...
    'k','linewidth',1);
hold on
plot3(src_loc(:,1),src_loc(:,2),src_loc(:,3),...
    'bo','linewidth',1,'markersize',5,'markerfacecolor','b');
plot3(tongue_loc(1),tongue_loc(2),tongue_loc(3),...
    'ro','linewidth',1,'markersize',5,'markerfacecolor','r');
axis equal
grid on
set(gca,'xtick',-5:5,'ytick',-5:5,'ztick',-5:5);
axis(axlim)

switch vh
    case 'h'
        view([-90 90]);
    case 'v'
        view([0 0]);
end
xlabel('X (cm)'); ylabel('Y (cm)'); zlabel('Z (cm)');
