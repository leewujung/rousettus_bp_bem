% 2016 03 14  Plot 3D BEM results with different mesh sizes

% folder_pre = 'C:\Users\Wu-Jung Lee\Dropbox\0_ANALYSIS\bp_processing';
folder_pre = 'F:\Dropbox\0_ANALYSIS\bp_processing';
fname_pre = 'flathead';

freq = 35*1e3;
% nn_all = [158, 341, 162, 210, 203];  % msh 0.2
% nn_all = [875, 1046, 1003, 906, 990];  % msh 0.2
% nn_all  = [389, 308, 310, 395, 825];  % msh 0.1
% nn_all = [3138, 3141, 3576, 3224, 3148];  % msh 0.1
% nn_all  = [3099, 3549, 1089, 3314, 2474];  % msh 0.05

mesh_all = [0.05,0.1,0.2];
nn_all = [3099, 3549, 1089, 3314, 2474;  % msh 0.05
    389, 308, 310, 395, 825;  % msh 0.1
    158, 341, 162, 210, 203];  % msh 0.2

fig = figure;
for iM=1:length(mesh_all)
    folder = ['20160314_bem3D_flathead',num2str(mesh_all(iM))];
    
    BEM_fname = sprintf('%s_freq%dkHz_bem_fp_mtx.mat',fname_pre,freq/1e3);
    BEM = load(fullfile(folder_pre,folder,BEM_fname),'bem_fp');  % load bem_fp for mesh detail
    
    for iN=1:length(nn_all)
        
        if iM==1
            fname = sprintf('%s_nn%05d_freq%dkHz.mat',fname_pre,nn_all(iM,iN),freq/1e3);
        else
            fname = sprintf('%s_nn%04d_freq%dkHz.mat',fname_pre,nn_all(iM,iN),freq/1e3);
        end
        FP = load(fullfile(folder_pre,folder,fname));
        
        pp = FP.bem_results.pfield_dB_norm;
        pp = pp+40+3;
        
        figure(fig)
        subplot(2,3,iN);
        polar(BEM.bem_fp.theta/180*pi,pp');
        hold on
        
        %     figure(fig2)
        %     plot(BEM.bem_fp.theta,pp);
        %     hold on
    end
end