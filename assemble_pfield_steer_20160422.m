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


% Individual src locations
% nn_all = [344, 337, 331, 325, 12];    % midline left
% nn_all = [251, 244, 238, 232, 13];    % midline right
% nn_all = [344, 337, 331, 325, 12,...  % midline left
%           251, 244, 238, 232, 13];    % midline right
% nn_all  = [1999, 2438, 1626, 2448, 2014, 2642, 1581, 2554, 1764, 2644,...   % left
%            1101,  520, 1398,  814, 1464,  543, 1125,  892, 1191,  442];     % right
nn_all  = [1999, 2438, 1626, 2448, 2014, 2642, 1581, 2554, 1764, 2644];     % left
% nn_all  = [1101,  520, 1398,  814, 1464,  543, 1125,  892, 1191,  442];     % right

freq = 35*1e3;
cvec = 0:-3:-30;
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm('eckert4');
mstruct = defaultm(mstruct);
num_freq = length(freq);
sm_len = 10;

% Tongue clicking location
tongue_loc_all = [...
    -0.02,-0.01,0;...
    -0.02,-0.008,0;...
    -0.02,-0.006,0;...
    -0.02,-0.004,0;...
    -0.02,-0.002,0;...
    -0.02,0.000,0;...
    -0.02,0.002,0;...
    -0.02,0.004,0;...
    -0.02,0.006,0;...
    -0.02,0.008,0;...
    -0.02,0.01,0];
num_tongue = size(tongue_loc_all,1);
colorset = jet(num_tongue);

% Load one sample to get info
indiv_src_file = sprintf('%s_nn%04d_freq%02dkHz.mat',...
    indiv_src_path,nn_all(1),freq/1e3);
load(fullfile(base_path,indiv_src_path,indiv_src_file));  % load in sample bem_results
load(fullfile(base_path,bem_results.src_path,bem_results.src_param_file));  % load in shape & param

% Main loop for assembling field
fig_bp_all = figure('units','normalized','outerposition',[0 0 1 1]);

fig_cntr_all = figure('position',[680,558,1023,420]);
subplot(121)  % left panel: src location
plot3(shape.nodesb(:,1),shape.nodesb(:,2),shape.nodesb(:,3),'.','markersize',5);
hold on
hsrc = plot3(shape.nodesb(nn_all,1),shape.nodesb(nn_all,2),shape.nodesb(nn_all,3),...
    'ko','linewidth',2);
axis equal
grid on
subplot(122)  % right panel: -3dB contour
axesm(map_proj);
axis off
framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-90 90]);
gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
colormap(jet(num_tongue))
colorbar('Ticks',linspace(0+1/num_tongue/2,1-1/num_tongue/2,num_tongue),...
    'TickLabels',{num2str(tongue_loc_all(:,2)*1e3)},'location','southoutside');
tightmap

c3db_xy = cell(num_tongue,1);
pp = cell(num_tongue,1);
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
    [x,y] = mfwdtran(mstruct,bem_results.theta,bem_results.phi);
    
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
    figure(fig_bp_all)
    % top plot: bp globe
    subplot(2,num_tongue,iT)
    axesm(map_proj);
    axis off
    contourf(-x,y,pp_plot,cvec(2:end));  % with contour line
    contour(-x,y,pp_plot,cvec(2:end),'fill','on');  % no contour line
    framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
    gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
    colorbar('Ticks',linspace(0+1/num_freq/2,1-1/num_freq/2,num_freq),...
        'TickLabels',{num2str(freq'/1e3)},'location','southoutside');
    tightmap
    colormap(parula(length(cvec)-1))
    caxis(cvec([end 1]))
    % bottom panel: src location
    subplot(2,num_tongue,iT+num_tongue)
    plot3(shape.nodesb(:,1),shape.nodesb(:,2),shape.nodesb(:,3),'.','markersize',5);
    hold on
    hsrc = plot3(shape.nodesb(nn_all,1),shape.nodesb(nn_all,2),shape.nodesb(nn_all,3),...
        'ko','linewidth',2);
    hton = plot3(tongue_loc(1),tongue_loc(2),tongue_loc(3),...
        'r*','linewidth',2,'markersize',10);
    axis equal
    grid on
    view([-90 90])
    
    % Plot -3dB contour figure
    figure(fig_cntr_all)
    subplot(121)
    % left panel: src location
    plot3(tongue_loc(1),tongue_loc(2),tongue_loc(3),...
        '*','linewidth',2,'markersize',10,'color',colorset(iT,:));
    hold on
    view([-90 90])
    % right panel: -3dB contour
    subplot(122)
    plot(-c3db_x,c3db_y,'linewidth',2,'color',colorset(iT,:));
    
    % Plot single tongue loc & bp
    fig_bp_single = figure('position',[680,558,1023,420]);
    subplot(121)  % left panel: src location
    plot3(shape.nodesb(:,1),shape.nodesb(:,2),shape.nodesb(:,3),'.','markersize',5);
    hold on
    hsrc = plot3(shape.nodesb(nn_all,1),shape.nodesb(nn_all,2),shape.nodesb(nn_all,3),...
        'ko','linewidth',2);
    hton = plot3(tongue_loc(1),tongue_loc(2),tongue_loc(3),...
        'r*','linewidth',2,'markersize',10);
    axis equal
    grid on
    view([-90 90])
    subplot(122);  % right panel: bp
    axesm(map_proj);
    axis off
%     contourf(-x,y,pp_plot,cvec(2:end));  % with contour line
    contour(-x,y,pp_plot,cvec(2:end),'fill','on');  % no contour line
    framem('fedgecolor',200*ones(1,3)/255,'flonlimit',[-180 180]);
    gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
    tightmap
    colormap(parula(length(cvec)-1))
    caxis(cvec([end 1]))
    colorbar('Ticks',sort(cvec),'Ticklabels',{num2str(sort(cvec)')},...
             'location','southoutside');

    title_txt = sprintf('tongue pos: (%d,%d,%d) mm',...
        tongue_loc(:,1)*1e3,tongue_loc(:,2)*1e3,tongue_loc(:,3)*1e3);
    suptitle(title_txt)

    if save_plot_opt
        save_bp_single_name = sprintf('%s_bp_tongue_%02d_%02d_%02d.png',...
            script_name,tongue_loc(:,1)*1e3,tongue_loc(:,2)*1e3,tongue_loc(:,3)*1e3);
        save_bp_single_name = regexprep(save_bp_single_name,'-','m');
        saveSameSize(fig_bp_single,'file',fullfile(save_path,save_bp_single_name),...
            'format','png','renderer','painters');
    end
    close(fig_bp_single)
end

if save_plot_opt
    save_cntr_name = sprintf('%s_mf_cntr.png',script_name);
    save_cntr_name = regexprep(save_cntr_name,'-','m');
    saveSameSize(fig_cntr_all,'file',fullfile(save_path,save_cntr_name),...
        'format','png','renderer','painters');
    
    save_bp_name = sprintf('%s_bp_all.png',script_name);
    save_bp_name = regexprep(save_bp_name,'-','m');
    saveSameSize(fig_bp_all,'file',fullfile(save_path,save_bp_name),...
        'format','png','renderer','painters');
end

%     close(fig_cntr_all)
%     close(fig_bp_all)

