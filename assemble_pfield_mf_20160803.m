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
az_plot_limit = 90;
num_freq = length(freq_all);
sm_len = 10;
% colorset = jet(num_freq);
cmap = '*Oranges';
colorset = brewermap(num_freq,cmap);

% Tongue clicking location
tongue_loc_all = [0.001,0.00,-0.00];

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
for iT=1%:num_tongue

    % Tonuge location
    tongue_loc = tongue_loc_all(iT,:);

    t_name = sprintf('%s_x%03d_y%03d_z%03d',...
        model_shape, tongue_loc(1)*1e3,tongue_loc(2)*1e3,tongue_loc(3)*1e3);
    suptitle_text = sprintf('Model %s, Tongue loc (%d,%d,%d) mm',...
        model_shape, tongue_loc(1)*1e3,tongue_loc(2)*1e3,tongue_loc(3)*1e3);

    % Prep figs
    fig_bp = figure('units','normalized','outerposition',[0 0 1 1]);
    
    fig_cntr = figure('position',[200,100,1400,420]);
    plot_model_shape_src(subplot(131),...  % shape, horizontal
        h_outline,shape.nodesb(nn_all,1:3),nan(1,3),h_axlim,'h');
    subplot(131)  % horizontal tonuge location
    plot3(tongue_loc(1)*1e2,tongue_loc(2)*1e2,tongue_loc(3)*1e2,...
        'x','linewidth',2,'markersize',8,'color','r');
    plot_model_shape_src(subplot(132),...  % shape, vertical
        v_outline,shape.nodesb(nn_all,1:3),nan(1,3),v_axlim,'v');
    subplot(132)  % vertical tongue location
    plot3(tongue_loc(1)*1e2,tongue_loc(2)*1e2,tongue_loc(3)*1e2,...
        'x','linewidth',2,'markersize',8,'color','r');
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
            indiv_src_file = sprintf('%s_nn%04d_freq%02dkHz.mat',...
                indiv_src_path,nn_all(iN),freq/1e3);
            load(fullfile(base_path,indiv_src_path,indiv_src_file));
            
            % phase delay
            gap_loc = shape.nodesb(nn_all(iN),1:3);
            src_gap_dist = sqrt((tongue_loc-gap_loc)*(tongue_loc-gap_loc)');
            phase_delay = exp(-1i*k*src_gap_dist);
            
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
    
        % Plot summary figure
        figure(fig_bp)
        plot_model_bp_globe(subplot(3,ceil((num_freq+2)/3),iF),...  % bp globe
            M,cvec,map_proj,az_plot_limit);
        title(sprintf('%02d kHz',freq/1e3));
    
        % Plot -3dB contour at current freq
        figure(fig_cntr)
        subplot(133)  % -3dB contour
        plot(-c3db_x,c3db_y,'linewidth',2,'color',colorset(iF,:));
        
    end  % loop through all freq

    figure(fig_bp)
    plot_model_shape_src(subplot(3,ceil((num_freq+2)/3),num_freq+1),...  % shape, horizontal
        h_outline,shape.nodesb(nn_all(:),1:3),tongue_loc,h_axlim,'h');
    subplot(3,ceil((num_freq+2)/3),num_freq+1);  % horizontal tonuge location
    plot3(tongue_loc(1)*1e2,tongue_loc(2)*1e2,tongue_loc(3)*1e2,...
        'x','linewidth',2,'markersize',8,'color','r');
    plot_model_shape_src(subplot(3,ceil((num_freq+2)/3),num_freq+2),...  % shape, vertical
        v_outline,shape.nodesb(nn_all(:),1:3),tongue_loc,v_axlim,'v');
    subplot(3,ceil((num_freq+2)/3),num_freq+2);  % vertical tongue location
    plot3(tongue_loc(1)*1e2,tongue_loc(2)*1e2,tongue_loc(3)*1e2,...
        'x','linewidth',2,'markersize',8,'color','r');
    
    suptitle(suptitle_text);
    save_bp_name = sprintf('%s_%s_bp_all',script_name,t_name);
    saveas(fig_bp,fullfile(save_path,[save_bp_name,'.fig']),'fig');
    saveSameSize(fig_bp,'file',fullfile(save_path,[save_bp_name,'.png']),...
        'format','png','renderer','painters');

    
    figure(fig_cntr)
    suptitle(suptitle_text);
    axesm(map_proj);
    axis off
    framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-1 1]*az_plot_limit);
    gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
%     colormap(jet(num_freq))
    colormap(brewermap(num_freq,cmap));
    colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
        'TickLabels',{num2str(freq_all'/1e3)},'location','southoutside');
    tightmap
    set(gcf,'color','w')
    
    save_cntr_name = sprintf('%s_%s_mf_cntr',script_name,t_name);
    saveas(fig_cntr,fullfile(save_path,[save_cntr_name,'.fig']),'fig');
    saveSameSize(fig_cntr,'file',fullfile(save_path,[save_cntr_name,'.png']),...
        'format','png','renderer','painters');

    close(fig_bp)
    close(fig_cntr)
end  % loop through tongue locations


