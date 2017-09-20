% 2016 03 14  Plot 3D BEM results with different mesh sizes

folder_pre = 'C:\Users\Wu-Jung Lee\Dropbox\0_ANALYSIS\bp_processing';

freq = 35*1e3;

% 2D BEM results
folder = '20160228_bem_10src';
fname_pre = 'head_src';
nn_all = [226, 254, 280, 306, 332];   % 2D BEM
fig = figure;
for iN=1:length(nn_all)
    fname = sprintf('%s_nn%d_freq_%dkHz.mat',fname_pre,nn_all(iN),freq/1e3);
    BEM2D = load(fullfile(folder_pre,folder,fname));
    
    pp = BEM2D.bem_output.pfield_dB_norm;
    pp = pp-(max(pp)-3);
    pp(pp<-40) = NaN;
    pp = pp+ 40;
   
    figure(fig)
    subplot(2,3,iN);
    h = polar(BEM2D.bem_output.phi,pp);
    set(h,'linewidth',2)
    hold on
end

% 3D BEM results
fname_pre = 'flathead';

% nn_all = [158, 341, 162, 210, 203];  % msh 0.2
% nn_all = [875, 1046, 1003, 906, 990];  % msh 0.2
% nn_all  = [389, 308, 310, 395, 825];  % msh 0.1
% nn_all = [3138, 3141, 3576, 3224, 3148];  % msh 0.1
% nn_all  = [3099, 3549, 1089, 3314, 2474];  % msh 0.05

mesh_all = [0.05,0.1,0.2];
nn_all = [3099, 3549, 1089, 3314, 2474;  % msh 0.05
          389, 308, 310, 395, 825;  % msh 0.1
          158, 341, 162, 210, 203];  % msh 0.2

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
        h = polar(BEM.bem_fp.theta/180*pi,pp');
        if iM==1
            set(h,'linewidth',2)
        end
        hold on
        title(['Gap ',num2str(iN)]);
    end
end




