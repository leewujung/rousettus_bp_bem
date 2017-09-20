function save_bem(A,B,CConst,save_path,save_bem_fname)

bem.A = A;
bem.B = B;
bem.CConst = CConst;

save(fullfile(save_path,save_bem_fname),'bem','-v7.3');