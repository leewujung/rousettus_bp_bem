function fp = load_fp(save_path,fp_fname)

A = load(fullfile(save_path,fp_fname));
fp = A.bem_fp;
