function plot_model_bp_globe_surf(h,M,cvec,map_proj,az_plot_range)


mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);

idx_nan = M.az>az_plot_range | M.az<-az_plot_range;
M.az(idx_nan) = NaN;
M.el(idx_nan) = NaN;
[M.x,M.y] = mfwdtran(mstruct,M.el,M.az);

M.az = -M.az;  % flip left/right because convention for azimuth is flipped in map projection
cgrey = 200*ones(1,3)/255;

axes(h)
axesm(map_proj);
framem('fedgecolor',cgrey,'flonlimit',[-1 1]*az_plot_range);
gridm('gcolor',cgrey,'glinestyle','-');
surfm(M.el,M.az,M.pp_plot);
tightmap
% colormap(parula(length(cvec)-1))
caxis(cvec([end 1]))
axis off
% colorbar('Ticks',sort(cvec),'Ticklabels',{num2str(sort(cvec)')},...
%     'location','southoutside');
