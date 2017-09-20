function plot_model_bp_globe(h,M,cvec,map_proj,az_plot_range)


mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);

idx_nan = M.az>az_plot_range | M.az<-az_plot_range;
M.az(idx_nan) = NaN;
M.el(idx_nan) = NaN;
[M.x,M.y] = mfwdtran(mstruct,M.el,M.az);

M.x = -M.x;  % flip left/right because convention for azimuth is flipped in map projection
cgrey = 200*ones(1,3)/255;

axes(h)
axesm(map_proj);
axis off
contourf(M.x,M.y,M.pp_plot,cvec(2:end),'w');  % with contour line
% contour(M.x,M.y,M.pp_plot,cvec(2:end),'fill','on');  % no contour line
gridm('gcolor',cgrey,'glinestyle','-');
framem('fedgecolor',cgrey,'flonlimit',[-1 1]*az_plot_range,...
    'flinewidth',1);
tightmap
colormap(parula(length(cvec)-1))
caxis(cvec([end 1]))
% colorbar('Ticks',sort(cvec),'Ticklabels',{num2str(sort(cvec)')},...
%     'location','southoutside');
