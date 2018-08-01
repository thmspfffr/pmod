%% FITTING
% pmod_wc_wholebrain_plot.m
Ies         = -4:0.1:-1;
Iis         = -5:0.1:-1;
Gg          = 0.62;
Gains       = -0.5:0.25:0.5;
% connectivity, AAL
v_conn =  1;
% simulations
v_sim = 1;
% dfa, aal, lcmv
v_dfa = 2;
% -------------------

Gs = 1;

% ---------
% LOAD EMPIRICAL DFA (LCMV)
% ---------
load(sprintf(['~/pupmod/proc/conn/' 'pupmod_src_dfa_v%d.mat'],v_dfa));
% ---------
% LOAD EMPIRICAL KURAMOTO PARAMETER
% ---------
% % subj x m x foi x cond
% load(sprintf(['~/pupmod/proc/conn/' 'pupmod_all_kuramoto_v%d.mat'],v));
% kura_emp_rest = mean(kura_std(:,1,1,1)./kura_mean(:,1,1,1));
% kura_emp_task = mean(kura_std(:,1,1,2)./kura_mean(:,1,1,2));
% ---------
% MATCH DFA WITH AAL ORDER USED FOR SIMULATIONS
% ---------
pars = [];
pars.grid = 'medium';
pars.N = 90;

dfa_emp_rest = nanmean(outp.dfa_all(:,:,1,1,1),2);
dfa_emp_task = nanmean(outp.dfa_all(:,:,1,1,2),2);

lambda_emp_rest = nanmean(outp.lambda_all(:,:,1,1,1),2);
lambda_emp_task = nanmean(outp.lambda_all(:,:,1,1,2),2);

clear dfa r dfa_r dist_fc fc_sim autocorr dfa_sim dfa_env_sim r*

load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v_conn));
mask = logical(tril(ones(90,90),-1));
mask = find(triu(ones(90))-eye(90));

fc_rest     =  squeeze(nanmean(cleandat(:,:,:,1,1,6),3));
fc_task     =  squeeze(nanmean(cleandat(:,:,:,1,2,6),3));
fc_rest_var =  std(nanmean(squeeze(nanmean(cleandat(:,:,:,1,1,6),3)))./max(nanmean(squeeze(nanmean(cleandat(:,:,:,1,1,6),3)))));

if ~exist(sprintf('~/pmod/proc/pmod_wholebrain_gainexcit_all_v%d.mat',v_sim))
  
  for ei = 1 : 2
    ei
    for igain = 1 : length(Gains)
      for iexc = 1 : length(Excit)
        %         igain
        load(sprintf('~/pmod/proc/pmod_wc_wholebrain_gainexcit_gain%d_excit%d_v%d.mat',igain,iexc,v_sim))
        
        if round(Ies(iies)*100)/100 == -2.8 && round( Iis(iiis)*100)/100 == -3.4
          disp('save stuff')
          FFCC = out.FC;
        end
        
        % Time scales
        outp.lambda(:,iies,iiis,iG,igain)      = mean(out.lambda,2);
        outp.lambda_env(:,iies,iiis,iG,igain)  = mean(out.lambda_env,2);
        % DFA
%         outp.dfa_sim(:,iies,iiis,iG,igain)     = squeeze(mean(out.dfa,1));
%         outp.dfa_env_sim(:,iies,iiis,iG,igain) = squeeze(mean(out.dfa_env,1));
        
        [~,peak_idx]=max(smooth(mean(out.PSD(out.f>3,:),2),20));
        outp.peakfreq(iies,iiis,iG,igain) = out.f(peak_idx+find(out.f<3,1,'last'));
        
        pars.dim = 2;
        outp.fc_sim_tmp = tp_match_aal(pars,mean(out.FC,3));
        outp.fc_sim_tmp = mean(out.FC,3);
        outp.fc_sim_var(:,iies,iiis,iG,igain) = std(nanmean(outp.fc_sim_tmp)./max(nanmean(outp.fc_sim_tmp)));
        
        [outp.r_rest_corr(iies,iiis,iG,igain), outp.p_rest_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_tmp(mask),fc_rest(mask));
        [outp.r_rest_corr_avg(iies,iiis,iG,igain), outp.p_rest_corr_avg(iies,iiis,iG,igain)]=corr(nanmean(outp.fc_sim_tmp)',nanmean(fc_rest)');
        
        outp.r_rest_corr_unc(iies,iiis,iG,igain) = dot(outp.fc_sim_tmp(mask),fc_rest(mask)) / sqrt(dot(outp.fc_sim_tmp(mask),outp.fc_sim_tmp(mask)) * dot(fc_rest(mask),fc_rest(mask)));
        %
        pars.dim = 1;
        
%         [outp.dfa_r_rest(iies,iiis,iG,igain), outp.dfa_p_rest(iies,iiis,iG,igain)] = corr(dfa_emp_rest(:),tp_match_aal(pars,repmat(outp.dfa_sim(:,iies,iiis,iG,igain),[1 90]),pars));
        %
        [outp.lambda_r_rest(iies,iiis,iG,igain), outp.lambda_p_rest(iies,iiis,iG,igain)] = corr(lambda_emp_rest,tp_match_aal(pars,repmat(outp.lambda(:,iies,iiis,iG,igain),[1 90]),pars));
        [outp.lambda_r_task(iies,iiis,iG,igain), outp.lambda_p_task(iies,iiis,iG,igain)] = corr(lambda_emp_task,tp_match_aal(pars,repmat(outp.lambda(:,iies,iiis,iG,igain),[1 90]),pars));
        
        outp.dist_fc_rest (iies,iiis,iG,igain)  = mean(outp.fc_sim_tmp(mask))-mean(fc_rest(mask));
        outp.dist_fc_task (iies,iiis,iG,igain)  = mean(outp.fc_sim_tmp(mask))-mean(fc_task(mask));
        
        outp.fc_sim_mean(iies,iiis,iG,igain) = mean(outp.fc_sim_tmp(mask));
        outp.fc_sim_all(:,:,iies,iiis,iG,igain) = outp.fc_sim_tmp;
        %         outp.fc_sim_env(iies,iiis,iG,igain) = mean(outp.fc_sim_env_tmp(mask));
        %
        outp.Ies(iies) = out.Ie;
        outp.Iis(iiis) = out.Ii;
        
        % KURAMOTO
        outp.kuramoto_mean (iies,iiis,iG,igain) = mean(out.KOPmean);
        outp.kuramoto_std (iies,iiis,iG,igain)  = mean(out.KOPsd);
        
        fclose all;
        %
      end
    end
  end
  
  
  % -----------------------------------
  % COMPUTE DEGREE IN MODEL
  % -----------------------------------
  
  vec = 1 : 90;
  
  for iG = 1 : Gs
    for igain = 1 : length(Gains)
      igain
      for ii = 1 : length(Iis)
        ii
        for ie = 1 : length(Ies)
          
          for ij = 1 : 90
            for jj = 1 : 90
              
              jjj = vec(vec~=ij);
              fc_tmp = outp.fc_sim_all(ij,jj,ie,ii,iG,igain);
              
              x_ref1 = mean(outp.fc_sim_all(ij,jjj,ie,ii,:,3),2);
              
              th(ij,jj,ie,ii,iG,igain) = fc_tmp>x_ref1;
            end
          end
          th_all(:,:,ie,ii,iG,igain) = th(:,:,ie,ii,iG,igain)+fliplr(rot90(th(:,:,ie,ii,iG,igain)>0,-1));
        end
      end
    end
    for igain = 1 : length(Gains)
      for ii = 1 : length(Iis)
        for ie = 1 : length(Ies)
          outp.degree(ie,ii,iG,igain) = nanmean(nanmean(th_all(:,:,ie,ii,iG,igain)));
        end
      end
    end
  end
  
  save(sprintf('~/pmod/proc/pmod_wholebrain_gainexcit_all_v%d.mat',v_sim),'outp')
  
else
  load(sprintf('~/pmod/proc/pmod_wholebrain_gainexcit_all_v%d.mat',v_sim))
end

error('!')

