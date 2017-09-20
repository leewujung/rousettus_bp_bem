function outline = get_head_shape_horizontal(nodes)
% Get shape outline from top view (looking -Z-axis)
% INPUT
%    nodes  head shape [Nx3]
% OUTPUT
%    outline  head outline looking down -Z-axis

k = boundary(nodes(:,1),nodes(:,2));
outline = nodes(k,:);

%% Below: using elevation angle
% [~,el,~] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));
% nn_ring = find(el>-1/180*pi & el<1/180*pi);  % for bullethead
% k = boundary(nodes(nn_ring,1),nodes(nn_ring,2));
% nn_ring = nn_ring(k);
% outline = nodes(nn_ring,:);
% nn_ring = [nn_ring;nn_ring(1,:)];
% outline = nodes(nn_ring,:);