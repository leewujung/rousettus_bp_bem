% 2015 03 14  3D BEM code directly from example
clear

% Set save path
usrn = getenv('username');
% save_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\20160315_bem3D_flathead0.1_surf'];
save_path = 'F:\Dropbox\0_ANALYSIS\bp_processing\20160315_bem3D_flathead0.1_surf';
save_fpre = 'flathead';

if ~exist(save_path,'dir')
    mkdir(save_path);
end

freq_all = (20:5:55)*1e3;
% freq_all = 35e3;


for iF=1:length(freq_all)
    
    save_fname = sprintf('%s_freq%dkHz_bem_fp_mtx.mat',save_fpre,freq_all(iF)/1e3);
    load(fullfile(save_path,save_fname))
    
    % Velocity
    nn_all  = [389, 308, 310, 395, 825, 3138, 3141, 3576, 3224, 3148];  % middle row
%     nn_all = [1017, 617, 903, 429, 709, 437, 479, 863, 824, 503];  % right side 2 rows
%     nn_all = [3198, 3654, 3537, 3766, 3566, 3201, 3711, 3282, 3610, 3677];  % left side 2 rows

    for iN=1:length(nn_all)
        
        nn = nn_all(iN);
        vp = zeros(shape.M,1);
        vp(nn) = param.u0;
        
        bem_results.src_index = nn;
        bem_results.vp = vp;
        
        % Pressure on surface
        pp = bem.A\(-bem.B*vp);
        
        bem_results.psurface = pp;
        
        % Pressure at field points
        ppFP = (bem_fp.Afp*pp+1j*param.k*param.rho*param.c*bem_fp.Bfp*vp)./bem_fp.Cfp;
        
        bem_results.pfield = ppFP;
        bem_results.pfield_dB = 20*log10(abs(ppFP));
        bem_results.pfield_dB_norm = bem_results.pfield_dB - max(bem_results.pfield_dB);
        
        bem_results.theta = bem_fp.theta;
        bem_results.phi = bem_fp.phi;
        bem_results.x = bem_fp.x;
        bem_results.y = bem_fp.y;
        bem_results.z = bem_fp.z;
        bem_results.xyzFP = bem_fp.xyzFP;

        save_fname = sprintf('%s_nn%04d_freq%dkHz.mat%',save_fpre,nn,freq_all(iF)/1e3);
        save(fullfile(save_path,save_fname),'bem_results');

%         figure(fig);
%         polar((-180:5:180)/180*pi,dd(10,:));
%         hold on
%         
%         figure;
%         axesm eckert4;
%         gridm('gcolor',190*ones(1,3)/255,'glinestyle','-');
%         surfm(theta,phi-90,dd);
%         caxis([0 30])
    end

end
