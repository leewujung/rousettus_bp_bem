% 2016 04 21  Calculate individual src pfield
% 2016 08 05  Revise for parallel computation and load src loc separately

% clear
usrn = getenv('username');

% Set up various paths
if strcmp(usrn,'Wu-Jung')
    base_path = 'F:\Dropbox\0_ANALYSIS\bp_bem_modeling';
else
    base_path = ['C:\Users\',usrn,'\Dropbox\0_ANALYSIS\bp_bem_modeling'];
end

freq_all = 35*1e3;
mtx_path = 'calc_bem_20160803_Ra224-0.6mm';
ss = strsplit(mtx_path,'_');
model_shape = ss{end};
bem_calc_date = ss{3};

[~,script_name,~] = fileparts(mfilename('fullpath'));
script_name = sprintf('%s_%s_%s',script_name,bem_calc_date,model_shape);
save_path = fullfile(base_path,script_name);
if ~exist(save_path,'dir')
    mkdir(save_path);
end

time_stamp = datestr(now,'yyyy-mm-dd, HH:MM:SS');  % timestamp

% Src locations path/file
src_loc_path = sprintf('pick_src_loc_%s',bem_calc_date);
src_loc_file = sprintf('%s_src_loc.mat',model_shape);
SRC = load(fullfile(base_path,src_loc_path,src_loc_file));
nn_all = [SRC.idx_right,SRC.idx_left];

% Load shape and gap locations
param_fname = sprintf('%s_freq%02dkHz_param.mat',mtx_path,35);
M = load(fullfile(base_path,mtx_path,param_fname));

parfor iF=1:length(freq_all)
    % Load pre-calculated bem stuff
    bem_fname = sprintf('%s_freq%02dkHz_bem_mtx.mat',mtx_path,freq_all(iF)/1e3);
    fp_fname = sprintf('%s_freq%02dkHz_fp_mtx.mat',mtx_path,freq_all(iF)/1e3);
    param_fname = sprintf('%s_freq%02dkHz_param.mat',mtx_path,freq_all(iF)/1e3);
    bem = load_bem(fullfile(base_path,mtx_path),bem_fname)
    bem_fp = load_fp(fullfile(base_path,mtx_path),fp_fname)
    [shape,param] = load_shape_param(fullfile(base_path,mtx_path),param_fname);
    
    for iN=1:length(nn_all)
        nn = nn_all(iN);
        vp = zeros(shape.M,1); vp(nn)=param.u0;
        
        % Pressure on surface
        pp = bem.A\(-bem.B*vp);
        
        % Pressure at field points
        ppFP = (bem_fp.Afp*pp+1j*param.k*param.rho*param.c*bem_fp.Bfp*vp)./bem_fp.Cfp;
        
        % Save results
        save_fname = sprintf('%s_nn%04d_freq%dkHz.mat%',script_name,nn,freq_all(iF)/1e3);
        %save(fullfile(save_path,save_fname),'bem_results');
        save_bem_results(save_path,save_fname,...
            time_stamp,src_loc_path,src_loc_file,nn,vp,...
            mtx_path,bem_fname,fp_fname,param_fname,...
            pp,ppFP,bem_fp);
    end
end

