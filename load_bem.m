function bem = load_bem(save_path,bem_fname)

A = load(fullfile(save_path,bem_fname));
bem = A.bem;
