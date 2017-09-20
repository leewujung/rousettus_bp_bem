function save_param(param,shape,k,save_path,save_param_fname)

param.k = k;

save(fullfile(save_path,save_param_fname),'param','shape');