% 2016 04 21  Assemble pfield from individual src
% 2016 05 07  Plot for NIFTI poster
% 2016 05 20  Plot for ASA talk

clear
usrn = getenv('username');

% Set up various paths
if strcmp(usrn,'Wu-Jung')
    base_path = 'F:\Dropbox\0_ANALYSIS\bp_bem_modeling';
    addpath('F:\Dropbox\0_CODE\beampattern_other_code');
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
else
    base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_bem_modeling'];
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\beampattern_other_code']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2']);
    addpath(['C:\Users\',usrn,'\Dropbox\0_CODE\MATLAB\saveSameSize']);
end

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

save_plot_opt = 1;

indiv_src_path = 'calc_indiv_src_pfield_20160421_bullethead';
ss = strsplit(indiv_src_path,'_');
model_shape = ss{end};

A.base_path = base_path;
A.indiv_src_path = indiv_src_path;
A.model_shape = model_shape;


% Load gap locations
mtx_path = 'calc_bem_20160421_bullethead';
param_file = sprintf('%s_freq%02dkHz_param.mat',mtx_path,35);
M = load(fullfile(base_path,mtx_path,param_file));
NN = M.shape.nodesb(:,1:3);
[az,el,r] = cart2sph(NN(:,1),NN(:,2),NN(:,3));
nn_rgape = find( ((el>7/180*pi & el<12/180*pi) |...
                 (el>-12/180*pi & el<-7/180*pi)) &...
                 az>-110/180*pi & az<0/180*pi);
nn_lgape = find( ((el>7/180*pi & el<12/180*pi) |...
                 (el>-12/180*pi & el<-7/180*pi)) &...
                 az>0/180*pi & az<110/180*pi);
nn_all = nn_lgape;

% Individual src locations
% nn_all = [344, 337, 331, 325, 12];    % midline left
% nn_all = [251, 244, 238, 232, 13];    % midline right
% nn_all = [344, 337, 331, 325, 12,...  % midline left
%           251, 244, 238, 232, 13];    % midline right
% nn_all  = [1999, 2438, 1626, 2448, 2014, 2642, 1581, 2554, 1764, 2644,...   % left
%            1101,  520, 1398,  814, 1464,  543, 1125,  892, 1191,  442];     % right
% nn_all  = [1999, 2438, 1626, 2448, 2014, 2642, 1581, 2554, 1764, 2644];     % left
% nn_all  = [1101,  520, 1398,  814, 1464,  543, 1125,  892, 1191,  442];     % right

freq_all = [25,35,45,55]*1e3;
cvec = 0:-3:-30;
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm('eckert4');
mstruct = defaultm(mstruct);
num_freq = length(freq_all);
colorset = jet(num_freq);
sm_len = 10;
az_plot_range = 120;

A.colorbar_vec = cvec;
A.freq_all = freq_all;
A.map.map_projection = map_proj;
A.map.mstruct = mstruct;
A.smooth_len = sm_len;
A.az_plot_range = az_plot_range;

% Tongue clicking location
tongue_loc_all = [-0.02,0.004,0.005];

% Load one sample to get info
indiv_src_file = sprintf('%s_nn%04d_freq%02dkHz.mat',...
    indiv_src_path,nn_all(1),freq_all(1)/1e3);
load(fullfile(base_path,indiv_src_path,indiv_src_file));  % load in sample bem_results

% Get corner extent
[xxx,yyy] = mfwdtran(mstruct,[0 90],[180 0]);
xxx = xxx(1);
yyy = yyy(2);

% Main loop for assembling field
c3db_xy = cell(num_freq,1);
pp = cell(num_freq,1);
for iT=1:size(tongue_loc_all,1)
    load(fullfile(base_path,bem_results.src_path,bem_results.src_param_file));  % load in shape & param
    
    tongue_loc = tongue_loc_all(iT,:);
    A.tongue_loc = tongue_loc;
    
    % Plot tongue and src location on model
    fig_bp = figure;
    set(fig_bp,'position',[100   100   900   600])
    if num_freq<=5
        subplot(num_freq,3,1:3:num_freq*3);
    else
        num_row_fig = round(num_freq/2);
        subplot(num_row_fig,3,1:3:num_row_fig*3);
    end
    plot3(shape.nodesb(:,1)*100,shape.nodesb(:,2)*100,shape.nodesb(:,3)*100,'.','markersize',5);
    hold on
    hsrc = plot3(shape.nodesb(nn_all,1)*100,shape.nodesb(nn_all,2)*100,shape.nodesb(nn_all,3)*100,...
        'ko','linewidth',2);
    hton = plot3(tongue_loc(1)*100,tongue_loc(2)*100,tongue_loc(3)*100,...
        'r*','linewidth',2,'markersize',10);
    xlabel('X (cm)')
    ylabel('Y (cm)')
    zlabel('Z (cm)')
    axis equal
    grid on
    title(model_shape)
%     tongue_legend = sprintf('tongue loc (%d,%d,%d) mm',...
%         tongue_loc(1)*1e3,tongue_loc(2)*1e3,tongue_loc(3)*1e3);
%     ll = legend([hsrc,hton],{'src loc',tongue_legend},'location','southoutside');
%     ll = legend([hsrc,hton],{'src loc','tongue loc'});
%     set(ll,'fontsize',12);
    view([-90 90]);
    
   
    for iF=1:num_freq
        
        k = 2*pi*freq_all(iF)/param.c;
        A.k = k;
        
        % Assemble total field
        pfield = nan(size(bem_results.pfield));
        for iN=1:length(nn_all)
            % load individual src pfield
            indiv_src_file = sprintf('%s_nn%04d_freq%02dkHz.mat',...
                indiv_src_path,nn_all(iN),freq_all(iF)/1e3);
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
        pp_plot(abs(bem_results.phi)>az_plot_range) = NaN;

        
        % Get -3dB main contour
        [~,vq_norm,azq,elq] = interp_bp(bem_results.phi/180*pi,bem_results.theta/180*pi,pp,'natural');
        azq = azq/pi*180;
        elq = elq/pi*180;
        idx = abs(azq)>90;
        vq_in = vq_norm;
        vq_in(idx) = NaN;
        [~,c_main_nan] = get_main_contour(vq_in,unique(azq),unique(elq),-3);  % get main contour at -3dB with NaN insert for break contour
        [c3db_x,c3db_y] = mfwdtran(mstruct,c_main_nan(:,2),c_main_nan(:,1));  % [az,el] to [x,y]
        
        % Plot summary figure
        figure(fig_bp)
        if num_freq<=5
            subplot(num_freq,3,(iF-1)*3+2:(iF-1)*3+3)
        else
            if iF<=num_row_fig
                subplot(num_row_fig,3,2+3*(iF-1));
            else
                subplot(num_row_fig,3,3+3*(iF-num_row_fig-1));
            end
        end
        
        [x,y] = mfwdtran(mstruct,bem_results.theta,bem_results.phi);
        axesm(map_proj);
%         contour(-x,y,pp_plot,cvec(2:end),'fill','on');
        contourf(-x,y,pp_plot,cvec(2:end),'fill','on','linecolor','w');
%         contourfm(elq/pi*180,azq/pi*180,vq_norm,contour_vec(2:cvec_min_idx),...
%             'fill','on','linecolor','w');  % don't plot 0 dB contour
        framem('fedgecolor',200*ones(1,3)/255,'flonlimit',az_plot_range*[-1 1]);
        gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
        axis off
        colormap(parula(length(cvec)-1))
        caxis(cvec([end 1]))
%         tightmap
        text(xxx-0.5,yyy-0.1,sprintf('%02dkHz',freq_all(iF)/1e3),'fontsize',12);
        
        fig_cntrf = figure;
        [x,y] = mfwdtran(mstruct,bem_results.theta,bem_results.phi);
        axesm(map_proj);
%         contour(-x,y,pp_plot,cvec(2:end),'fill','on');
        contourf(-x,y,pp_plot,cvec(2:end),'fill','on','linecolor','w');
%         contourfm(elq/pi*180,azq/pi*180,vq_norm,contour_vec(2:cvec_min_idx),...
%             'fill','on','linecolor','w');  % don't plot 0 dB contour
        framem('fedgecolor',200*ones(1,3)/255,'flonlimit',az_plot_range*[-1 1]);
        gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
        axis off
        colormap(parula(length(cvec)-1))
        caxis(cvec([end 1]))
%         tightmap
        text(xxx-0.5,yyy-0.1,sprintf('%02dkHz',freq_all(iF)/1e3),'fontsize',12);
        tightmap
        
    end
    
    title_txt = sprintf('tongue pos: (%d,%d,%d) mm',...
        tongue_loc(:,1)*1e3,tongue_loc(:,2)*1e3,tongue_loc(:,3)*1e3);

    figure(fig_bp)
    suptitle(title_txt);
    
    
    if save_plot_opt
        save_all_name = sprintf('%s_summary_tongue_%02d_%02d_%02d',...
            script_name,tongue_loc(:,1)*1e3,tongue_loc(:,2)*1e3,tongue_loc(:,3)*1e3);
        save_all_name = regexprep(save_all_name,'-','m');
        saveas(fig_bp,fullfile(save_path,[save_all_name,'.fig']),'fig');
        saveSameSize(fig_bp,'file',fullfile(save_path,save_all_name),...
            'format','png','renderer','painters');
        
        save_contrf_name = sprintf('%s_contourf_tongue_%02d_%02d_%02d',...
            script_name,tongue_loc(:,1)*1e3,tongue_loc(:,2)*1e3,tongue_loc(:,3)*1e3);
        save_contrf_name = regexprep(save_contrf_name,'-','m');
        saveas(fig_cntrf,fullfile(save_path,[save_contrf_name,'.fig']),'fig');
        saveSameSize(fig_cntrf,'file',fullfile(save_path,save_contrf_name),...
            'format','png','renderer','painters');
    end
%     close(fig_bp)
end


