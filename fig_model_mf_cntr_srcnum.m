% 2016 08 12  Plot model bp multi-freq contorus
%             Revised from assemble_pfield_mf_20160803
% 2017 10 05  Add code to find beam centers and plot them

clear
% Set up various paths
if isunix
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/rbfinterp_v1.2');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/saveSameSize');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/EllipseDirectFit');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/brewermap');
    save_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp/bp_bem_modeling';
    model_base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp/bp_bem_modeling/';
else
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    addpath('F:\Dropbox\0_CODE\MATLAB\brewermap');
    save_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp\bp_bem_modeling';
    model_base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp\bp_bem_modeling';
end

% Individual pfield path
% indiv_src_path = 'calc_indiv_src_pfield_20160829_Ra-colony-noear-0.5mm';
% indiv_src_path = 'calc_indiv_src_pfield_20160829_Ra-colony-0.7mm';
% indiv_src_path = 'calc_indiv_src_pfield_20160820_bullethead-sc-colony-0.5mm';
% indiv_src_path = 'calc_indiv_src_pfield_20160817_Ra-colony-0.5mm';
% indiv_src_path = 'calc_indiv_src_pfield_20160811_Ra224-0.6mm';
% indiv_src_path = 'calc_indiv_src_pfield_20160812_bullethead-sc-0.5mm';
indiv_src_path = 'calc_indiv_src_pfield_20160917_Ra-colony-rotear-0.5mm';
ss = strsplit(indiv_src_path,'_');
model_shape = ss{end};
bem_calc_date = ss{end-1};
pfield_size = [73,37];
lr_opt = 'right';

% Save path
[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(save_base_path,sprintf('%s_%s_%s',script_name,bem_calc_date,model_shape));
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Src locations path/file
src_incl_str = '3456';
src_incl = [3,4,5,6];
% src_incl2 = [1,3,4,6];
src_loc_path = sprintf('pick_src_loc_%s_%s',bem_calc_date,model_shape);
src_loc_file = sprintf('%s_src_loc_newleft34.mat',model_shape);
% src_loc_file = sprintf('%s_src_loc.mat',model_shape);
% src_loc_file = sprintf('%s_src_loc_many.mat',model_shape);
SRC = load(fullfile(model_base_path,src_loc_path,src_loc_file));
if strcmp(lr_opt,'left')
    nn_all = SRC.idx_left;
else
    nn_all = SRC.idx_right;
end

% Set other params
% freq_str = '20-24.9-25-25.1-30-35-40';
% freq_all = [20,24.9,25,25.1,30,35,40]*1e3;
freq_str = '25-5-55';
freq_all = [25:5:55]*1e3;
% freq_all(1) = 26;
% freq_all = [24.2:0.2:24.8]*1e3;
% freq_all = 35*1e3;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tongue clicking location
% tongue_loc = [0.028,0.000,-0.003];  % Used for bullethead-sc-colony-0.5mm
tongue_loc = [0.029,0.000,-0.005];  % Used for Ra-colony-rotear-0.5mm
% tongue_loc = [0.027,0.000,-0.005];  % Used for Ra-colony-noear-0.5mm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_fname = sprintf('%s_%s_%s_x%05.1f_y%05.1f_z%05.1f_%s_%s_f%skHz',...
    script_name,bem_calc_date,model_shape,...
    tongue_loc(1,1)*1e3,tongue_loc(1,2)*1e3,tongue_loc(1,3)*1e3,...
    src_incl_str,lr_opt,freq_str);
title_text = sprintf('%s %s, tongue (%5.1f,%5.1f,%5.1f), freq %s kHz',...
    bem_calc_date,model_shape,...
    tongue_loc(1,1)*1e3,tongue_loc(1,2)*1e3,tongue_loc(1,3)*1e3,...
    freq_str);

% Load one sample to get info
if rem(freq_all(1)/1e3,1)==0
    indiv_src_file = sprintf('%s_nn%04d_freq%dkHz.mat',...
        indiv_src_path,nn_all(1),freq_all(1)/1e3);
else
    indiv_src_file = sprintf('%s_nn%04d_freq%2.1fkHz.mat',...
        indiv_src_path,nn_all(1),freq_all(1)/1e3);
end
load(fullfile(model_base_path,indiv_src_path,indiv_src_file));  % load in sample bem_results
load(fullfile(model_base_path,bem_results.src_path,bem_results.src_param_file));  % load in shape & param

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
    h_outline,shape.nodesb(nn_all,1:3),tongue_loc,h_axlim,'h');
% subplot(131)  % horizontal tonuge location
% plot3(tongue_loc(1)*1e2,tongue_loc(2)*1e2,tongue_loc(3)*1e2,...
%     'x','linewidth',2,'markersize',8,'color','r');
plot_model_shape_src(subplot(132),...  % shape, vertical
    v_outline,shape.nodesb(nn_all,1:3),tongue_loc,v_axlim,'v');
% subplot(132)  % vertical tongue location
% plot3(tongue_loc(1)*1e2,tongue_loc(2)*1e2,tongue_loc(3)*1e2,...
%     'x','linewidth',2,'markersize',8,'color','r');
subplot(133)  % -3dB contour
axesm(map_proj);
axis off
hold on

for iF=1:num_freq
    
    freq = freq_all(iF);
    k = 2*pi*freq/param.c;
    
    % Assemble total field
   pfield = nan(size(bem_results.pfield));
    for iN=1:length(nn_all)
        % load individual src pfield
        if rem(freq/1e3,1)==0
            indiv_src_file = sprintf('%s_nn%04d_freq%dkHz.mat',...
                indiv_src_path,nn_all(iN),freq/1e3);
        else
            indiv_src_file = sprintf('%s_nn%04d_freq%2.1fkHz.mat',...
                indiv_src_path,nn_all(iN),freq/1e3);
        end
        load(fullfile(model_base_path,indiv_src_path,indiv_src_file));
        
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
    
    % Find beam center
    idxnotnan = ~isnan(M.pp_plot);
    [~,M.vq_norm,M.azq,M.elq] = ...
        interp_bp(M.az(idxnotnan)/180*pi,M.el(idxnotnan)/180*pi,M.pp_plot(idxnotnan),'rbf');
    M.azq = M.azq/pi*180;
    M.elq = M.elq/pi*180;

    % --- max beam energy location
    xx = M.vq_norm(:);
    [~,mm_idx] = max(xx);
    M.max.el = M.elq(mm_idx);
    M.max.az = M.azq(mm_idx);
    % --- averaged location of normalized beam energy >-1
    xx(isnan(xx)) = -inf;
    [~,sort_idx] = sort(xx,'descend');
    ii = xx(sort_idx)>-1;
    M.top.el = mean(M.elq(sort_idx(ii)));
    M.top.az = mean(M.azq(sort_idx(ii)));
    % --- fit ellipse
    % ------- M.el = bem_results.theta
    % ------- M.az = bem_results.phi
    [raw,rot_max,rot_elpctr,rot_elpctr_tilt] = ...
        shift_rotate_bp_composite(M.az,M.el,M.pp,map_proj,0.005);
    [el_ectr,az_ectr] = minvtran(mstruct,rot_max.E.x0,rot_max.E.y0);  % inverse map projection
    [M.ectr.el,M.ectr.az] = rotatem(el_ectr,az_ectr,...
                                      [M.max.el,M.max.az],...
                                      'inverse','degrees');

    % Plot summary figure
    figure(fig_bp)
    plot_model_bp_globe(subplot(3,ceil((num_freq+2)/3),iF),...  % bp globe
        M,cvec,map_proj,az_plot_limit);
    hold on
    plotm(M.max.el,-M.max.az,'rx','markersize',8,'linewidth',2);
    plotm(M.top.el,-M.top.az,'r^','markersize',8,'linewidth',2);
    plotm(M.ectr.el,-M.ectr.az,'ro','markersize',8,'linewidth',2);
    if rem(freq/1e3,1)==0
        title(sprintf('%02d kHz',freq/1e3));
    else
        title(sprintf('%2.1f kHz',freq/1e3));
    end
    
    % Plot -3dB contour at current freq
    figure(fig_cntr)
    subplot(133)  % -3dB contour
    plot(-c3db_x,c3db_y,'linewidth',2,'color',colorset(iF,:));
    
end  % loop through all freq

figure(fig_bp)
plot_model_shape_src(subplot(3,ceil((num_freq+2)/3),num_freq+1),...  % shape, horizontal
    h_outline,shape.nodesb(nn_all(:),1:3),tongue_loc,h_axlim,'h');
plot_model_shape_src(subplot(3,ceil((num_freq+2)/3),num_freq+2),...  % shape, vertical
    v_outline,shape.nodesb(nn_all(:),1:3),tongue_loc,v_axlim,'v');

suptitle(title_text);
% saveas(fig_bp,fullfile(save_path,[save_fname,'_all.fig']),'fig');
saveSameSize(fig_bp,'file',fullfile(save_path,[save_fname,'_all.png']),...
    'format','png','renderer','painters');
epswrite(fullfile(save_path,[save_fname,'_all.eps']));


figure(fig_cntr)
suptitle(title_text);
axesm(map_proj);
axis off
framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-1 1]*180);
gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
colormap(jet(num_freq))
% colormap(brewermap(num_freq,cmap));
colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
    'TickLabels',{num2str(freq_all'/1e3)},'location','southoutside');
tightmap
set(gcf,'color','w')

% saveas(fig_cntr,fullfile(save_path,[save_fname,'_cntr.fig']),'fig');
saveSameSize(fig_cntr,'file',fullfile(save_path,[save_fname,'_cntr.png']),...
    'format','png','renderer','painters');
epswrite(fullfile(save_path,[save_fname,'_cntr.eps']));

% close(fig_bp)
% close(fig_cntr)



