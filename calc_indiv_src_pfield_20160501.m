% 2016 04 21  Calculate individual src pfield

clear
usrn = getenv('username');

% Set up various paths
if strcmp(usrn,'Wu-Jung')
    base_path = 'F:\Dropbox\0_ANALYSIS\bp_bem_modeling';
else
    base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_bem_modeling'];
end

freq_all = [20:5:60]*1e3;
mtx_path = 'calc_bem_20160501_bullethead';
ss = strsplit(mtx_path,'_');
model_shape = ss{end};

[~,script_name,~] = fileparts(mfilename('fullpath'));
script_name = [script_name,'_',model_shape];
save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

% Load gap locations
mtx_path = 'calc_bem_20160501_bullethead';
param_file = sprintf('%s_freq%02dkHz_param.mat',mtx_path,35);
M = load(fullfile(base_path,mtx_path,param_file));
NN = M.shape.nodesb(:,1:3);
nn_all = [358,666,464,732,704,4309,1479,237,962,520,...        % right
         4273,2147,2461,2512,369,1602,332,2597,2384,2401];  % left

for iF=1:length(freq_all)
    % Load pre-calculated bem stuff
    bem_file = sprintf('%s_freq%02dkHz_bem_mtx.mat',mtx_path,freq_all(iF)/1e3);
    fp_file = sprintf('%s_freq%02dkHz_fp_mtx.mat',mtx_path,freq_all(iF)/1e3);
    param_file = sprintf('%s_freq%02dkHz_param.mat',mtx_path,freq_all(iF)/1e3);
    load(fullfile(base_path,mtx_path,bem_file));
    load(fullfile(base_path,mtx_path,fp_file));
    load(fullfile(base_path,mtx_path,param_file));
    bem_results.src_path = mtx_path;
    bem_results.src_bem_file = bem_file;
    bem_results.src_fp_file = fp_file;
    bem_results.src_param_file = param_file;
    
    for iN=1:length(nn_all)
        nn = nn_all(iN);
        vp = zeros(shape.M,1); vp(nn)=param.u0;
        
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

        % Save results
        save_fname = sprintf('%s_nn%04d_freq%dkHz.mat%',script_name,nn,freq_all(iF)/1e3);
        save(fullfile(save_path,save_fname),'bem_results');
    end
end

