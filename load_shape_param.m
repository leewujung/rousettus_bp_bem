function [shape,param] = load_shape_param(save_path,param_fname)

A = load(fullfile(save_path,param_fname));
shape = A.shape;
param = A.param;