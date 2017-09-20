function save_bem_results(save_path,save_fname,...
            time_stamp,src_loc_path,src_loc_file,nn,vp,...
            mtx_path,bem_fname,fp_fname,param_fname,...
            pp,ppFP,bem_fp)

bem_results.time_stamp = time_stamp;     
bem_results.src_loc_path = src_loc_path;
bem_results.src_loc_file = src_loc_file;

bem_results.src_path = mtx_path;
bem_results.src_bem_file = bem_fname;
bem_results.src_fp_file = fp_fname;
bem_results.src_param_file = param_fname;

bem_results.src_index = nn;
bem_results.vp = vp;
bem_results.psurface = pp;

bem_results.pfield = ppFP;
bem_results.pfield_dB = 20*log10(abs(ppFP));
bem_results.pfield_dB_norm = bem_results.pfield_dB - max(bem_results.pfield_dB);

bem_results.theta = bem_fp.theta;
bem_results.phi = bem_fp.phi;
bem_results.x = bem_fp.x;
bem_results.y = bem_fp.y;
bem_results.z = bem_fp.z;
bem_results.xyzFP = bem_fp.xyzFP;

save(fullfile(save_path,save_fname),'bem_results');
