% 2016 04 21  Assemble pfield from individual src
% 2016 07 27  Clean up code
% 2016 08 16  Plot for paper; steer bp horizontally
% 2017 09 20  Plot center of the beam on all figures
%             (both max energy point and center of best-fitting ellipse)
%             Also return those locations for quantification

clear
usrn = getenv('username');
% Set up various paths
if isunix
    base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp/bp_bem_modeling';
    addpath('~/code_git/rousettus_bp');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/rbfinterp_v1.2');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/saveSameSize');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/brewermap');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/mtit');
else
    base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp\bp_bem_modeling';
    addpath('F:\Dropbox\0_CODE\beampattern_other_code');
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    addpath('F:\Dropbox\0_CODE\MATLAB\brewermap');
    addpath('F:\Dropbox\0_CODE\MATLAB\mtit');
end

% Individual pfield path
indiv_src_path = 'calc_indiv_src_pfield_20160917_Ra-colony-rotear-0.5mm';
% indiv_src_path = 'calc_indiv_src_pfield_20160820_bullethead-sc-colony-0.5mm';
% indiv_src_path = 'calc_indiv_src_pfield_20160817_Ra-colony-0.5mm';
% indiv_src_path = 'calc_indiv_src_pfield_20160811_Ra224-0.6mm';
% indiv_src_path = 'calc_indiv_src_pfield_20160812_bullethead-sc-0.5mm';
ss = strsplit(indiv_src_path,'_');
model_shape = ss{end};
bem_calc_date = ss{end-1};
pfield_size = [73,37];
lr_opt = 'left';

% Save path
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(base_path,sprintf('%s_%s_%s',script_name,bem_calc_date,model_shape));
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Src locations path/file
src_incl_str = '2345';
src_incl = [2,3,4,5];
src_loc_path = sprintf('pick_src_loc_%s_%s',bem_calc_date,model_shape);
src_loc_file = sprintf('%s_src_loc_newleft34.mat',model_shape);
SRC = load(fullfile(base_path,src_loc_path,src_loc_file));
if strcmp(lr_opt,'left')
    nn_all = SRC.idx_left;
else
    nn_all = SRC.idx_right;
end

% Set other params
freq = 35*1e3;
cvec = 0:-3:-30;
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);
az_plot_limit = 90;
num_freq = length(freq);
sm_len = 10;

% Tongue clicking location
% tongue_loc_all = ...  % Used for bullethead-sc-colony-0.5mm
%    [0.029,0.000,-0.003;...
%     0.028,0.000,-0.003;...
%     0.027,0.000,-0.003;...
%     0.026,0.000,-0.003;...
%     0.025,0.000,-0.003];
% tongue_loc_all = ...  % Used for Ra-colony-rotear-0.5mm
%    [0.029,0.002,-0.005;...
%     0.028,0.002,-0.005;...
%     0.027,0.002,-0.005;...
%     0.026,0.002,-0.005;...
%     0.025,0.002,-0.005];
tongue_loc_all = ...  % Used for Ra-colony-rotear-0.5mm
   [0.029,0.000,-0.005;...
    0.0275,0.000,-0.005;...
    0.026,0.000,-0.005;...
    0.0245,0.000,-0.005;...
    0.023,0.000,-0.005];
tongue_loc_str = '029t023';

save_fname = sprintf('%s_%s_%s_x%s_y%05.1f_z%05.1f_%s_%ddeg_%s',...
    script_name,bem_calc_date,model_shape,...
    tongue_loc_str,tongue_loc_all(1,2)*1e3,tongue_loc_all(1,3)*1e3,...
    src_incl_str,az_plot_limit,lr_opt);
title_text = sprintf('%s %s, tongue (%s,%5.1f,%5.1f)',...
    bem_calc_date,model_shape,...
    regexprep(tongue_loc_str,'t','-'),tongue_loc_all(1,2)*1e3,tongue_loc_all(1,3)*1e3);

% Tongue clicking location
num_tongue = size(tongue_loc_all,1);
cmap = '*Reds';
cmap_ctr = '*Blues';
% colorset = brewermap(num_tongue,cmap);
colorset = brewermap(num_tongue+1,cmap);
colorset_ctr = brewermap(num_tongue+3,cmap_ctr);

% Load one sample to get info
indiv_src_file = sprintf('%s_nn%04d_freq%02dkHz.mat',...
    indiv_src_path,nn_all(1),freq/1e3);
load(fullfile(base_path,indiv_src_path,indiv_src_file));  % load in sample bem_results
load(fullfile(base_path,bem_results.src_path,bem_results.src_param_file));  % load in shape & param

% Sort and select nn_all
if ~isempty(src_incl)
    nn_all_x = shape.nodesb(nn_all,1);
    [~,idx_x_sort] = sort(nn_all_x,'descend');
    nn_all = nn_all(idx_x_sort);  % sort according x locations
    nn_all = nn_all(src_incl);  % include only selected sources
end

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

%%

% Prep figs
fig_bp = figure('units','normalized','outerposition',[0 0 1 1]);

fig_cntr = figure('position',[200,100,1400,420]);
plot_model_shape_src(subplot(131),...  % shape, horizontal
    h_outline,shape.nodesb(nn_all,1:3),nan(1,3),h_axlim,'h');
plot_model_shape_src(subplot(132),...  % shape, vertical
    v_outline,shape.nodesb(nn_all,1:3),nan(1,3),v_axlim,'v');
subplot(133)  % -3dB contour
axesm(map_proj);
axis off
hold on


%% Loop through tongue position
for iT=1:num_tongue
    
    tongue_loc = tongue_loc_all(iT,:);
    k = 2*pi*freq/param.c;
    
    % Assemble total field
    pfield = nan(size(bem_results.pfield));
    for iN=1:length(nn_all)
        % load individual src pfield
        indiv_src_file = sprintf('%s_nn%04d_freq%02dkHz.mat',...
            indiv_src_path,nn_all(iN),freq/1e3);
        load(fullfile(base_path,indiv_src_path,indiv_src_file));
        
        % phase delay
        gap_loc = shape.nodesb(nn_all(iN),1:3);
        src_gap_dist = sqrt((tongue_loc-gap_loc)*(tongue_loc-gap_loc)');
        phase_delay = exp(1i*k*src_gap_dist);
        
        pfield(:,iN) = bem_results.pfield*phase_delay;
    end
    
    pfield_tot = sum(pfield,2);
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
    
    % Get -3dB main contour
    [~,vq_norm,azq,elq] = interp_bp(bem_results.phi/180*pi,bem_results.theta/180*pi,pp,'natural');
    azq = azq/pi*180;
    elq = elq/pi*180;
    idx = abs(azq)>90;
    vq_in = vq_norm;
    vq_in(idx) = NaN;
    [~,c_main_nan] = get_main_contour(vq_in,unique(azq),unique(elq),-3);  % get main contour at -3dB with NaN insert for break contour
    [c3db_x,c3db_y] = mfwdtran(mstruct,c_main_nan(:,2),c_main_nan(:,1));  % [az,el] to [x,y]
    c3db_x = smooth(c3db_x,sm_len);
    c3db_y = smooth(c3db_y,sm_len);
    
    % Get beam center location
    % --- max beam energy location
    [~,mmidx] = max(vq_in(:));
    azq_max_loc = azq(mmidx);
    elq_max_loc = elq(mmidx);    
    % --- center of best-fitting ellipse
    [raw,rot_max,rot_elpctr,rot_elpctr_tilt] = ...
        shift_rotate_bp_composite(bem_results.phi,bem_results.theta,pp,map_proj,0.005);
    [el_ectr,az_ectr] = minvtran(mstruct,rot_max.E.x0,rot_max.E.y0);  % inverse map projection
    [el_ectr_r,az_ectr_r] = rotatem(el_ectr,az_ectr,...
                                    [elq_max_loc,azq_max_loc],...
                                    'inverse','degrees');
    
    % Plot summary figure
    figure(fig_bp)
    plot_model_bp_globe(subplot(3,num_tongue,iT),...  % bp globe
        M,cvec,map_proj,az_plot_limit);
    plot_model_shape_src(subplot(3,num_tongue,iT+num_tongue),...  % shape, horizontal
        h_outline,shape.nodesb(nn_all,1:3),tongue_loc,h_axlim,'h');
    plot_model_shape_src(subplot(3,num_tongue,iT+num_tongue*2),...  % shape, vertical
        v_outline,shape.nodesb(nn_all,1:3),tongue_loc,v_axlim,'v');
    subplot(3,num_tongue,iT);
    hold on
    % flip left/right because convention for azimuth is flipped in map projection
    plotm(el_ectr_r,-az_ectr_r,'ro','markersize',8,'linewidth',2);
    plotm(elq_max_loc,-azq_max_loc,'rx','markersize',8,'linewidth',2);
    
    % Plot -3dB contour at current freq
    figure(fig_cntr)
    subplot(131)  % horizontal
    plot3(tongue_loc(1)*1e2,tongue_loc(2)*1e2,tongue_loc(3)*1e2,...
        '.','linewidth',2,'markersize',20,'color',colorset(iT,:));
    subplot(132)  % vertical
    plot3(tongue_loc(1)*1e2,tongue_loc(2)*1e2,tongue_loc(3)*1e2,...
        '.','linewidth',2,'markersize',20,'color',colorset(iT,:));
    subplot(133)  % -3dB contour
    plot(-c3db_x,c3db_y,'linewidth',2,'color',colorset(iT,:));
    hold on
    % flip left/right because convention for azimuth is flipped in map projection
    plotm(el_ectr_r,-az_ectr_r,'.','color',colorset(iT,:),'markersize',30,'linewidth',2);
    plotm(elq_max_loc,-azq_max_loc,'x','color',colorset_ctr(iT,:),'markersize',8,'linewidth',2);
    
end

% Save -3dB contour figure
figure(fig_cntr)
subplot(133)
axesm(map_proj);
axis off
gridm('gcolor',200*ones(1,3)/255,'glinestyle','-');
framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-1 1]*az_plot_limit,...
    'flinewidth',1);
% colormap(brewermap(num_tongue,cmap));
% colorbar('Ticks',linspace(0+1/num_tongue/2,1-1/num_tongue/2,num_tongue),...
%     'TickLabels',{num2str(tongue_loc_all(:,1)*1e3)},'location','southoutside');
colormap(brewermap(num_tongue+1,cmap));
colorbar('Ticks',linspace(0+1/(num_tongue+1)/2,1-1/(num_tongue+1)/2,(num_tongue+1)),...
    'TickLabels',{num2str([tongue_loc_all(:,1);0]*1e3)},'location','southoutside');
tightmap
mtit(fig_cntr,title_text);  % call this to avoid chaning fig position

saveSameSize(fig_cntr,'file',fullfile(save_path,[save_fname,'_cntr.png']),...
    'format','png','renderer','painters');
epswrite(fullfile(save_path,[save_fname,'_cntr.eps']));

% Save all bp figure
figure(fig_bp)
suptitle(title_text);
saveSameSize(fig_bp,'file',fullfile(save_path,[save_fname,'_all.png']),...
    'format','png','renderer','painters');
epswrite(fullfile(save_path,[save_fname,'_all.eps']));


