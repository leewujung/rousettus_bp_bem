function outline = get_head_shape_vertical(nodes)
% Get shape outline from top view (looking -Z-axis)
% INPUT
%    nodes  head shape [Nx3]
% OUTPUT
%    outline  head outline looking down -Z-axis

k = boundary(nodes(:,1),nodes(:,3));
outline = nodes(k,:);

%% Below: use azimuth angle
% [az,~,~] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));
% nn_ring = find((az>-1/180*pi & az<1/180*pi) |...
%                (az>179/180*pi | az<-179/180*pi) );  % for bullethead
% k = boundary(nodes(nn_ring,1),nodes(nn_ring,3));
% nn_ring = nn_ring(k);
% nn_ring = [nn_ring;nn_ring(1,:)];
% outline = nodes(nn_ring,:);