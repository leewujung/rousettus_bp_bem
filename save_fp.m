function save_fp(save_path,save_fp_fname,...
                 time_stamp,param_fname,theta,phi,x,y,z,xyzFP,...
                 Afp,Bfp,Cfp)

bem_fp.code_start_time = time_stamp;
bem_fp.param_fname = param_fname;  % param file used
bem_fp.theta = theta;
bem_fp.phi = phi;
bem_fp.x = x;
bem_fp.y = y;
bem_fp.z = z;
bem_fp.xyzFP = xyzFP;
bem_fp.Afp = Afp;
bem_fp.Bfp = Bfp;
bem_fp.Cfp = Cfp;
save(fullfile(save_path,save_fp_fname),'bem_fp','-v7.3');