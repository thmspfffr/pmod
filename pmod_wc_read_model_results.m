function outp = pmod_wc_read_model_results(para)




for iies = 1 : length(Ies)
    for iiis = 1 : length(Iis)
      for iG = 1:length(Gg)
        for igain = 1:length(Gains)
          %         igain
%           fn = 
          load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',iies,iiis,iG,igain,para_v))

          outp.peakfreq(iies,iiis,iG,igain) = out.peakfreq;

          pars.dim = 2;
          
          outp.fc_sim_tmp = out.FC;
          outp.fc_sim_env_tmp = out.FC_env;
                    
          [outp.r_rest_corr(iies,iiis,iG,igain), outp.p_rest_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_tmp(mask),fc_rest(mask));
          [outp.r_env_rest_corr(iies,iiis,iG,igain), outp.p_env_rest_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_rest(mask));
%           [outp.r_env_task_corr(iies,iiis,iG,igain), outp.p_env_task_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_task(mask));
          
          [outp.r_rest_corr_avg(iies,iiis,iG,igain), outp.p_rest_corr_avg(iies,iiis,iG,igain)]=corr(nanmean(outp.fc_sim_tmp)',nanmean(fc_rest)');
          [outp.r_env_rest_indiv_corr(:,iies,iiis,iG,igain), outp.p_env_rest_indiv_corr(:,iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_rest_indiv(mask,:));
          %
          pars.dim = 1;
          [~,~,outp.kds(iies,iiis,iG,igain)] = kstest2(fc_rest(mask),outp.fc_sim_env_tmp(mask));
   
          outp.dist_fc_rest (iies,iiis,iG,igain)  = mean(outp.fc_sim_tmp(mask))-mean(fc_rest(mask));
          %         outp.dist_fc_task (iies,iiis,iG,igain)  = mean(outp.fc_sim_tmp(mask))-mean(fc_task(mask));
          
          outp.fc_sim_mean(iies,iiis,iG,igain)    = mean(outp.fc_sim_tmp(mask));
          outp.fc_sim_env_mean(iies,iiis,iG,igain) = mean(outp.fc_sim_env_tmp(mask));
          outp.Ies(iies) = out.Ie;
          outp.Iis(iiis) = out.Ii;
          
          fclose all;
          %
        end
      end
    end
  end
  
  
  %
  save(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_all_v%d.mat',v_sim),'outp')
%   % DELETE OUTOUT FILES
%   for iies = 1 : length(Ies)
%     iies
%     for iiis = 1 : length(Iis)
%       %     iiis
%       for iG = 1:length(Gg)
%         for igain = 1:length(Gains)
%           
%           fn =  sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d',iies,iiis,iG,igain,v_sim);
%           delete(sprintf('%s_processing.txt',fn))
%           delete(sprintf('%s.mat',fn))
%           
%         end
%       end
%     end
%   end
% else
%   load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_all_v%d.mat',v_sim))
% end

error('!')