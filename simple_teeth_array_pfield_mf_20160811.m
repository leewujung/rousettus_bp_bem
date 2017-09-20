% 2016 04 21  Assemble pfield from individual src
% 2016 07 27  Clean up code

clear
usrn = getenv('username');
% Set up various paths
if strcmp(usrn,'Wu-Jung')
    base_path = 'F:\Dropbox\0_ANALYSIS\bp_bem_modeling';
    addpath('F:\Dropbox\0_CODE\beampattern_other_code');
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    addpath('F:\Dropbox\0_CODE\MATLAB\brewermap');
else
    base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_bem_modeling'];
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\beampattern_other_code']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\brewermap']);
end

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end


indiv_src_path = 'calc_indiv_src_pfield_20160803_Ra224-0.5mm';
ss = strsplit(indiv_src_path,'_');
model_shape = ss{end};
bem_calc_date = ss{end-1};
pfield_size = [73,37];


% Src locations path/file
src_loc_path = sprintf('pick_src_loc_%s_%s',bem_calc_date,model_shape);
src_loc_file = sprintf('%s_src_loc.mat',model_shape);
SRC = load(fullfile(base_path,src_loc_path,src_loc_file));
nn_all = SRC.idx_left;

% Set other params
freq_all = [20:5:60]*1e3;
cvec = 0:-3:-30;
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);
az_plot_limit = 180;
num_freq = length(freq_all);
sm_len = 10;
colorset = jet(num_freq);
% cmap = '*Oranges';
% colorset = brewermap(num_freq,cmap);

% Tongue clicking location
tongue_loc = [0.008,-0.000,0];

% Load one sample to get info
indiv_src_file = sprintf('%s_nn%04d_freq%02dkHz.mat',...
    indiv_src_path,nn_all(1),freq_all(1)/1e3);
load(fullfile(base_path,indiv_src_path,indiv_src_file));  % load in sample bem_results
load(fullfile(base_path,bem_results.src_path,bem_results.src_param_file));  % load in shape & param

% Get shape outline
pad = 0.002;
v_outline = get_head_shape_vertical(shape.nodesb(:,1:3));
v_axlim = [min(v_outline(:,1))-pad,...
    max(v_outline(:,1))+pad,...
    -0.1,0.1,...
    min(v_outline(:,3))-pad,...
    max(v_outline(:,3))+pad];

h_outline = get_head_shape_horizontal(shape.nodesb(:,1:3));
h_axlim = [min(h_outline(:,1))-pad,...
    max(h_outline(:,1))+pad,...
    min(h_outline(:,2))-pad,...
    max(h_outline(:,2))+pad,...
    -0.1,0.1];


%% Loop through tongue position

x_pos = [3, 5.5, 9, 11, 15,...
    2, 4.5, 7.5, 10, 14]'*1e-3*1.2;
z_pos = [0, 0, 0, 0, 0,...
    -3.5, -4, -5, -6, -6]'*1e-3;
x_pos = x_pos-mean(x_pos);
z_pos = z_pos-mean(z_pos);
y_pos = zeros(size(x_pos));

rot_z_angle = -20/180*pi;
Rz = [cos(rot_z_angle), -sin(rot_z_angle), 0;...
      sin(rot_z_angle),  cos(rot_z_angle), 0;...
      0, 0, 1];
xyz_rot = (Rz*[x_pos,y_pos,z_pos]')';
x_pos = xyz_rot(:,1);
y_pos = xyz_rot(:,2)+0.005;
z_pos = xyz_rot(:,3);

[x,y,z] = sph2cart(bem_results.phi/180*pi,bem_results.theta/180*pi,1);
% x_pos = shape.nodesb(nn_all,1);
% y_pos = shape.nodesb(nn_all,2);
% z_pos = shape.nodesb(nn_all,3);
rn = sqrt((repmat(x',length(x_pos),1)-repmat(x_pos,1,length(x))).^2 +...
    (repmat(y',length(y_pos),1)-repmat(y_pos,1,length(y))).^2 +...
    (repmat(z',length(z_pos),1)-repmat(z_pos,1,length(z))).^2);
r = sqrt((x'-x_pos(1)).^2+(y'-y_pos(1)).^2+(z'-z_pos(1)).^2);
r_diff = rn-repmat(r,size(rn,1),1);
fac1 = 1./rn;

steer_d = sqrt((tongue_loc(1)-x_pos).^2+(tongue_loc(2)-y_pos).^2+(tongue_loc(3)-z_pos).^2);
steer_d = steer_d-min(steer_d);


% Tonuge location
t_name = sprintf('%s_x%03d_y%03d_z%03d',...
    model_shape, tongue_loc(1)*1e3,tongue_loc(2)*1e3,tongue_loc(3)*1e3);
suptitle_text = sprintf('Model %s, Tongue loc (%d,%d,%d) mm',...
    model_shape, tongue_loc(1)*1e3,tongue_loc(2)*1e3,tongue_loc(3)*1e3);

% Prep figs
fig_bp = figure('units','normalized','outerposition',[0 0 1 1]);

for iF=1:num_freq
    
    freq = freq_all(iF);
    k = 2*pi*freq/param.c;
    
    % Assemble total field
    fac2 = exp( 1j*(k*(r_diff+repmat(steer_d(:),1,size(r_diff,2)))) );
    %fac2 = exp(1j*k*r_diff);
    %fac3 = repmat(w_pos(:),1,size(r_diff,2));
    %b = sum(fac1.*fac2.*fac3,1);
    pfield_tot = sum(fac1.*fac2,1);
    
    pfield_tot_dB = 20*log10(abs(pfield_tot));
    pfield_tot_dB = pfield_tot_dB-max(max(pfield_tot_dB));
    pp = reshape(pfield_tot_dB,size(bem_results.phi));
    pp_plot = pp;
    pp_plot(pp_plot<cvec(end)) = cvec(end);
    
    pp = reshape(pp,pfield_size);
    pp_plot = reshape(pp_plot,pfield_size);
    bem_results.phi = reshape(bem_results.phi,pfield_size);
    bem_results.theta = reshape(bem_results.theta,pfield_size);
    
    M.el = bem_results.theta;
    M.az = bem_results.phi;
    M.pp = pp;  % original assembled pfield
    M.pp_plot = pp_plot;  % assembled pfield with values smaller than cvec(end) removed
    
    % Plot summary figure
    figure(fig_bp)
    plot_model_bp_globe(subplot(3,ceil((num_freq+2)/3),iF),...  % bp globe
        M,cvec,map_proj,az_plot_limit);
    title(sprintf('%02d kHz',freq/1e3));
    
%     % Plot -3dB contour at current freq
%     figure(fig_cntr)
%     subplot(133)  % -3dB contour
%     plot(-c3db_x,c3db_y,'linewidth',2,'color',colorset(iF,:));
    
end  % loop through all freq

figure(fig_bp)
plot_model_shape_src(subplot(3,ceil((num_freq+2)/3),num_freq+1),...  % shape, horizontal
    h_outline,[x_pos,y_pos,z_pos],tongue_loc,h_axlim,'h');
subplot(3,ceil((num_freq+2)/3),num_freq+1);  % horizontal tonuge location
plot3(tongue_loc(1)*1e2,tongue_loc(2)*1e2,tongue_loc(3)*1e2,...
    'x','linewidth',2,'markersize',8,'color','r');
plot_model_shape_src(subplot(3,ceil((num_freq+2)/3),num_freq+2),...  % shape, vertical
    v_outline,[x_pos,y_pos,z_pos],tongue_loc,v_axlim,'v');
subplot(3,ceil((num_freq+2)/3),num_freq+2);  % vertical tongue location
plot3(tongue_loc(1)*1e2,tongue_loc(2)*1e2,tongue_loc(3)*1e2,...
    'x','linewidth',2,'markersize',8,'color','r');

suptitle(suptitle_text);
save_bp_name = sprintf('%s_%s_bp_all',script_name,t_name);
saveas(fig_bp,fullfile(save_path,[save_bp_name,'.fig']),'fig');
saveSameSize(fig_bp,'file',fullfile(save_path,[save_bp_name,'.png']),...
    'format','png','renderer','painters');


