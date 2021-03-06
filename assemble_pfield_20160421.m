% 2016 04 21  Assemble pfield from individual src

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

% Individual src locations
% nn_all = [344, 337, 331, 325, 12];    % midline left
% nn_all = [251, 244, 238, 232, 13];    % midline right
% nn_all = [344, 337, 331, 325, 12,...  % midline left
%           251, 244, 238, 232, 13];    % midline right
% nn_all  = [1999, 2438, 1626, 2448, 2014, 2642, 1581, 2554, 1764, 2644,...   % left
%            1101,  520, 1398,  814, 1464,  543, 1125,  892, 1191,  442];     % right
nn_all  = [1999, 2438, 1626, 2448, 2014, 2642, 1581, 2554, 1764, 2644];     % left
% nn_all  = [1101,  520, 1398,  814, 1464,  543, 1125,  892, 1191,  442];     % right

freq_all = [20:5:60]*1e3;
cvec = 0:-3:-30;
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm('eckert4');
mstruct = defaultm(mstruct);
num_freq = length(freq_all);
colorset = jet(num_freq);
sm_len = 10;

A.colorbar_vec = cvec;
A.freq_all = freq_all;
A.map.map_projection = map_proj;
A.map.mstruct = mstruct;
A.smooth_len = sm_len;

% Tongue clicking location
tongue_loc_all = [...
    -0.02,0.000,0;...
    -0.02,0.002,0;...
    -0.02,0.004,0;...
    -0.02,0.006,0;...
    -0.02,0.008,0;...
    -0.02,0.01,0;...
    -0.02,-0.002,0;...
    -0.02,-0.004,0;...
    -0.02,-0.006,0;...
    -0.02,-0.008,0;...
    -0.02,-0.01,0];
% tongue_loc_all = [...
%     -0.02,0.005,0.000;...
%     -0.02,0.005,0.002;...
%     -0.02,0.005,0.004;...
%     -0.02,0.005,0.006;...
%     -0.02,0.005,0.008;...
%     -0.02,0.005,0.01;...
%     -0.02,0.005,-0.01;...
%     -0.02,0.005,-0.008;...
%     -0.02,0.005,-0.006;...
%     -0.02,0.005,-0.004;...
%     -0.02,0.005,-0.002];
% tongue_loc_all = 0.005;
% tongue_loc_all = [-0.02*ones(length(tongue_loc_all),1),...
%     tongue_loc_all',...
%     0*ones(length(tongue_loc_all),1)];

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
    fig_all = figure;
    set(fig_all,'position',[100   100   900   600])
    if num_freq<=5
        subplot(num_freq,3,1:3:num_freq*3);
    else
        num_row_fig = round(num_freq/2);
        subplot(num_row_fig,3,1:3:num_row_fig*3);
    end
    plot3(shape.nodesb(:,1),shape.nodesb(:,2),shape.nodesb(:,3),'.','markersize',5);
    hold on
    hsrc = plot3(shape.nodesb(nn_all,1),shape.nodesb(nn_all,2),shape.nodesb(nn_all,3),...
        'ko','linewidth',2);
    hton = plot3(tongue_loc(1),tongue_loc(2),tongue_loc(3),...
        'r*','linewidth',2,'markersize',10);
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    axis equal
    grid on
    title(model_shape)
%     tongue_legend = sprintf('tongue loc (%d,%d,%d) mm',...
%         tongue_loc(1)*1e3,tongue_loc(2)*1e3,tongue_loc(3)*1e3);
%     ll = legend([hsrc,hton],{'src loc',tongue_legend},'location','southoutside');
%     ll = legend([hsrc,hton],{'src loc','tongue loc'});
%     set(ll,'fontsize',12);
    view([-90 90]);
    
    % Plot multi-frequency contour
    fig_mf = figure;
    axesm(map_proj);
    axis off
    framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-90 90]);
    gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
    colormap(jet(num_freq))
    colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
        'TickLabels',{num2str(freq_all'/1e3)},'location','southoutside');
    tightmap
    
    
    
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
        %         pfield_tot_dB(pfield_tot_dB<-40) = NaN;
        %         pfield_tot_dB = pfield_tot_dB+40;
        pp = reshape(pfield_tot_dB,size(bem_results.phi));
        pp_plot = pp;
        pp_plot(pp_plot<cvec(end)) = cvec(end);
        
        A.pfield_indiv{iF} = pfield;     % individual pfield
        A.pfield_total{iF} = pfield_tot;
        A.pfield_total_dB{iF} = pfield_tot_dB;
        A.pfield_total_dB_2D{iF} = pp;
        A.pfield_total_dB_2D_forplot{iF} = pp_plot;
        A.lat{iF} = bem_results.theta;   % elevation
        A.lon{iF} = bem_results.phi;     % azimuth, need to flip +/-for plotting
        % to get the left/right direction right
        
        % Get -3dB main contour
        [~,vq_norm,azq,elq] = interp_bp(bem_results.phi/180*pi,bem_results.theta/180*pi,pp,'natural');
        azq = azq/pi*180;
        elq = elq/pi*180;
        idx = abs(azq)>90;
        vq_in = vq_norm;
        vq_in(idx) = NaN;
        [~,c_main_nan] = get_main_contour(vq_in,unique(azq),unique(elq),-3);  % get main contour at -3dB with NaN insert for break contour
        [c3db_x,c3db_y] = mfwdtran(mstruct,c_main_nan(:,2),c_main_nan(:,1));  % [az,el] to [x,y]
        
        A.contour_3dB_azel{iF} = c_main_nan;
        A.contour_3dB_xy{iF} = [c3db_x,c3db_y];
        
        % Plot summary figure
        figure(fig_all)
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
        contour(-x,y,pp_plot,cvec(2:end),'fill','on');
        framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
        gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
        axis off
        colormap(parula(length(cvec)-1))
        caxis(cvec([end 1]))
        tightmap
        text(xxx-0.5,yyy-0.1,sprintf('%02dkHz',freq_all(iF)/1e3),'fontsize',12);
        %         colorbar('southoutside','ticks',fliplr(cvec));
        
        %         figure(fig_all)
        %         subplot(num_freq,3,(iF-1)*3+2:(iF-1)*3+3)
        %         contourf(lon/pi*180-90,lat/pi*180,pp,cvec(2:end));
        %         grid on
        %         xlim([-180 180])
        %         set(gca,'xtick',-180:30:180,'ytick',-90:30:90);
        %         colormap(parula(length(cvec)-1))
        %         colorbar
        %         caxis(cvec([end 1]))
        
        % Plot multi-frequency contours
        c3db_x = smooth(c3db_x,sm_len);
        c3db_y = smooth(c3db_y,sm_len);
        figure(fig_mf)
        plot(-c3db_x,c3db_y,'linewidth',2,'color',colorset(iF,:));
    end
    
    title_txt = sprintf('tongue pos: (%d,%d,%d) mm',...
        tongue_loc(:,1)*1e3,tongue_loc(:,2)*1e3,tongue_loc(:,3)*1e3);

    figure(fig_all)
    suptitle(title_txt);
    
    figure(fig_mf);
%     framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
%     gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
    colormap(jet(num_freq))
    colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
        'TickLabels',{num2str(freq_all'/1e3)},'location','southoutside');
    tightmap
    title(title_txt);
    
    if save_plot_opt
        save_mf_name = sprintf('%s_mf_cntr_tongue_%02d_%02d_%02d.png',...
            script_name,tongue_loc(:,1)*1e3,tongue_loc(:,2)*1e3,tongue_loc(:,3)*1e3);
        save_mf_name = regexprep(save_mf_name,'-','m');
        saveSameSize(fig_mf,'file',fullfile(save_path,save_mf_name),...
            'format','png','renderer','painters');
        
        save_all_name = sprintf('%s_summary_tongue_%02d_%02d_%02d.png',...
            script_name,tongue_loc(:,1)*1e3,tongue_loc(:,2)*1e3,tongue_loc(:,3)*1e3);
        save_all_name = regexprep(save_all_name,'-','m');
        saveSameSize(fig_all,'file',fullfile(save_path,save_all_name),...
            'format','png','renderer','painters');
    end
    close(fig_all)
    close(fig_mf)
end


