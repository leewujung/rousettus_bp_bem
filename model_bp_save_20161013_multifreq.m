% 2016 04 21  Assemble pfield from individual src
% 2016 07 27  Clean up code
% 2016 08 04  Save model bp for projection sampling test
% 2017 10 05  Do multi-freq modeling version of bullethead model

clear
usrn = getenv('username');
% Set up various paths
if isunix
    base_path = '~/internal_2tb/Dropbox/Z_wjlee/projects/rousettus_bp/bp_bem_modeling';
    addpath('~/code_git/rousettus_bp');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/rbfinterp_v1.2');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/saveSameSize');
    addpath('~/internal_2tb/Dropbox/0_CODE/MATLAB/brewermap');
else
    base_path = 'F:\Dropbox\Z_wjlee\projects\rousettus_bp\bp_bem_modeling';
    addpath('F:\Dropbox\0_CODE\beampattern_other_code');
    addpath('F:\Dropbox\0_CODE\MATLAB\rbfinterp_v1.2');
    addpath('F:\Dropbox\0_CODE\MATLAB\saveSameSize');
    addpath('F:\Dropbox\0_CODE\MATLAB\brewermap');
end

[~,script_name,~] = fileparts(mfilename('fullpath'));
save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

indiv_src_path = 'calc_indiv_src_pfield_20160820_bullethead-sc-colony-0.5mm';
ss = strsplit(indiv_src_path,'_');
model_shape = ss{end};
bem_calc_date = ss{end-1};
pfield_size = [73,37];

S.base_path = base_path;
S.indiv_src_path = indiv_src_path;
S.model_shape = model_shape;
S.bem_calc_date = bem_calc_date;


% Src locations path/file
src_loc_path = sprintf('pick_src_loc_%s_%s',bem_calc_date,model_shape);
src_loc_file = sprintf('%s_src_loc.mat',model_shape);
SRC = load(fullfile(base_path,src_loc_path,src_loc_file));
nn_all = SRC.idx_left([4,3,1,2]);  % source 2-5 counting from the front

% Other params
freq_all = (20:5:60)*1e3;
%freq = 35*1e3;
cvec = 0:-3:-30;
map_proj = 'eckert4';   % use eckert4 because orthographic has very limited azimuth span
mstruct = defaultm(map_proj);
mstruct = defaultm(mstruct);
az_plot_limit = 180;
num_freq = length(freq_all);
sm_len = 10;

S.nn_all = nn_all;


% Tongue clicking location
tongue_loc_all = [0.028,0.000,-0.003];
num_tongue = size(tongue_loc_all,1);

for iF=1:num_freq
    freq = freq_all(iF);
    S.freq = freq;


    % Load one sample to get info
    indiv_src_file = sprintf('%s_nn%04d_freq%02dkHz.mat',...
                             indiv_src_path,nn_all(1),freq/1e3);
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
    for iT=1:num_tongue
        
        tongue_loc = tongue_loc_all(iT,:);
        k = 2*pi*freq/param.c;
        
        S.tongue_loc = tongue_loc;
        save_fname = sprintf('%s_%s_%02dkHz_x%03d_y%03d_z%03d',...
                             script_name,model_shape,freq/1e3,...
                             tongue_loc(1)*1e3,tongue_loc(2)*1e3,tongue_loc(3)*1e3);
        fig_title = sprintf('%s %02dkHz, tongue at (%d,%d,%d) mm',...
                            model_shape,freq/1e3,tongue_loc(1)*1e3,tongue_loc(2)*1e3,tongue_loc(3)*1e3);

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
        
        S.el = bem_results.theta;
        S.az = bem_results.phi;
        S.pp = pp;  % original assembled pfield
        S.pp_plot = pp_plot;  % assembled pfield with values smaller than cvec(end) removed

        
        % Plot bp and head shape
        fig_bp = figure('position',[200,100,1400,420]);
        plot_model_shape_src(subplot(131),...  % shape, horizontal
                             h_outline,shape.nodesb(nn_all,1:3),nan(1,3),h_axlim,'h');
        plot_model_shape_src(subplot(132),...  % shape, vertical
                             v_outline,shape.nodesb(nn_all,1:3),nan(1,3),v_axlim,'v');
        plot_model_bp_globe(subplot(133),...  % bp globe
                            S,cvec,map_proj,az_plot_limit);
        %     suptitle(fig_title);
        
        % Plot tongue location
        subplot(131)  % horizontal
        hold on
        plot3(tongue_loc(1)*1e2,tongue_loc(2)*1e2,tongue_loc(3)*1e2,...
              'rx','linewidth',2,'markersize',8);
        hold off
        subplot(132)  % vertical
        hold on
        plot3(tongue_loc(1)*1e2,tongue_loc(2)*1e2,tongue_loc(3)*1e2,...
              'rx','linewidth',2,'markersize',8);
        hold off
        
        % Save figure and mat file
        saveSameSize(fig_bp,'file',fullfile(save_path,[save_fname,'.png']),...
                     'format','png','renderer','painters');
        save(fullfile(save_path,[save_fname,'.mat']),'-struct','S');

        %     close(fig_bp)

    end


end  % loop through all freq