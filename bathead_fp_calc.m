% 2016 03 08  Use existing BEM mtx to calculate pressure at field points

clear

% Set path
usrn = getenv('username');
bem_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\20160308_bem3D_bullethead_0.1'];
save_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_processing\20160308_bem3D_bullethead_0.1_fp'];
if ~exist(save_path,'dir')
    mkdir(save_path);
end

fname_pre = 'head_src';

freq_all = (20:5:45)*1e3;
nn_all = [1206, 1223, 1079, 1260, 992];

for iF = 1:length(freq_all)
    fprintf('Processing %dkHz field points\n',freq_all(iF)/1e3);
    
    % Load pre-calculated BEM stuff
    bem_file = sprintf('%s_%dkHz_bem_mtx.mat',fname_pre,freq_all(iF)/1e3);
    D = load(fullfile(bem_path,bem_file));
    
    for iN = 1:length(nn_all)
        % Assign source
        vp = zeros(D.shape.M,1);
        vp(nn_all(iN)) = D.param.u0;
        
        % Solve for pressure on surface
        pp = D.bem.A\(-D.bem.B*vp);
        
        % Solve for pressure at field points
        ppFP = (D.bem_fp.Afp*pp + 1j*D.param.k*D.param.rho*D.param.c*D.bem_fp.Bfp*vp)./D.bem_fp.Cfp;
        ppFP = (         Afp*pp + 1j*        k*        rho*        c*         Bfp*vp)./         Cfp;
        % Prep for saving
        S.bem_path = bem_path;
        S.bem_file = bem_file;
        S.shape = D.shape;
        S.param = D.param;
        S.theta_fp = D.bem_fp.theta;
        S.phi_fp = D.bem_fp.phi;
        S.v_src = vp;
        S.p_surf = pp;
        S.p_fp = ppFP;
        
        % Save results
        save_fname = sprintf('%s_nn%04d_freq_%dkHz.mat%',fname_pre,nn_all(iN),freq_all(iF)/1e3);
        save(fullfile(save_path,save_fname),'-struct','S');
    end
    
    clear D
end





