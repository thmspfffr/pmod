%% FITTING
% pmod_wc_wholebrain_plot.m

%-------------------------------------------------------------------------
% VERSION 2: COMPUTE PEAK FREQ
%-------------------------------------------------------------------------
% v_sim           = 1;
% Ies         = -10:0.5:10;
% Iis         = -10:0.5:10;
% Gg          = 0:0.1:1;
% Gains       = 0;
% nTrials     = 1;
% tmax        = 10000; % in units of tauE
% wins        = [3 50];
%-------------------------------------------------------------------------
% VERSION 2: 
%-------------------------------------------------------------------------
v_sim           = 2;
Ies         = -5:0.5:7.5;
Iis         = -10:0.5:2.5;
Gg          = 0.7;
Gains       = 0:0.25:0.25;
nTrials     = 1;
tmax        = 10000; % in units of tauE
wins        = [3 50];
%-------------------------------------------------------------------------

% connectivity, AAL
v_conn =  1;
% simulations
% dfa, aal, lcmv
v_dfa = 2;
% -------------------

Gs = 1;
v = v_sim
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
cnt = 0;
fc_rest     =  squeeze(nanmean(cleandat(:,:,:,1,1,6),3));
fc_task     =  squeeze(nanmean(cleandat(:,:,:,1,2,6),3));
fc_rest_var =  std(nanmean(squeeze(nanmean(cleandat(:,:,:,1,1,6),3)))./max(nanmean(squeeze(nanmean(cleandat(:,:,:,1,1,6),3)))));
%%
if ~exist(sprintf('~/pmod/proc/pmod_wholebrain_woolrich_all_v%d.mat',v_sim))

for iies = 1: length(Ies)
  iies
  for iiis = 1: length(Iis)
    for iG = 1 : length(Gg)
      for igain = 1 : length(Gains)
        %
        %         igain
        try
          load(sprintf('~/pmod/proc/pmod_wc_woolrich_rest_Ie%d_Ii%d_G%d_gain%d_v%d.mat',iies,iiis,iG,igain,v))
          
          %         if round(Ies(iies)*100)/100 == -2.8 && round( Iis(iiis)*100)/100 == -3.4
          %           disp('save stuff')
          %           FFCC = out.FC;
          %         end
          
          % Time scales
          outp.lambda(:,iies,iiis,iG,igain)      = mean(out.lambda,2);
          outp.lambda_env(:,iies,iiis,iG,igain)  = mean(out.lambda_env,2);
          % DFA
          %         outp.dfa_sim(:,iies,iiis,iG,igain)     = squeeze(mean(out.dfa,1));
          %         outp.dfa_env_sim(:,iies,iiis,iG,igain) = squeeze(mean(out.dfa_env,1));
          
          
          outp.peakfreq(iies,iiis,iG,igain) = out.peakfreq;
          
          pars.dim = 2;
          outp.fc_sim_tmp = tp_match_aal(pars,mean(out.FC,3));
          outp.fc_sim_var(:,iies,iiis,iG,igain) = std(nanmean(outp.fc_sim_tmp)./max(nanmean(outp.fc_sim_tmp)));
          
          outp.fc_sim_env = tp_match_aal(pars,mean(out.FC_env,3));
          
          [outp.r_rest_corr(iies,iiis,iG,igain), outp.p_rest_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_tmp(mask),fc_rest(mask));
          [outp.r_rest_corr_avg(iies,iiis,iG,igain), outp.p_rest_corr_avg(iies,iiis,iG,igain)]=corr(nanmean(outp.fc_sim_tmp)',nanmean(fc_rest)');
          
          [outp.r_env_rest_corr(iies,iiis,i,igainG), outp.p_env_rest_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_env(mask),fc_rest(mask));
          
          outp.r_rest_corr_unc(iies,iiis,iG,igain) = dot(outp.fc_sim_tmp(mask),fc_rest(mask)) / sqrt(dot(outp.fc_sim_tmp(mask),outp.fc_sim_tmp(mask)) * dot(fc_rest(mask),fc_rest(mask)));
          %
          pars.dim = 1;
          
          %         [outp.dfa_r_rest(iies,iiis,iG,igain), outp.dfa_p_rest(iies,iiis,iG,igain)] = corr(dfa_emp_rest(:),tp_match_aal(pars,repmat(outp.dfa_sim(:,iies,iiis,iG,igain),[1 90]),pars));
          %
          [outp.lambda_r_rest(iies,iiis,iG,igain), outp.lambda_p_rest(iies,iiis,iG,igain)] = corr(lambda_emp_rest,tp_match_aal(pars,repmat(outp.lambda(:,iies,iiis,iG),[1 90]),pars));
          [outp.lambda_r_task(iies,iiis,iG,igain), outp.lambda_p_task(iies,iiis,iG,igain)] = corr(lambda_emp_task,tp_match_aal(pars,repmat(outp.lambda(:,iies,iiis,iG),[1 90]),pars));
          
          outp.dist_fc_rest (iies,iiis,iG,igain)  = mean(outp.fc_sim_tmp(mask))-mean(fc_rest(mask));
          outp.dist_fc_task (iies,iiis,iG,igain)  = mean(outp.fc_sim_tmp(mask))-mean(fc_task(mask));
          
          outp.fc_sim_mean(iies,iiis,iG,igain) = mean(outp.fc_sim_tmp(mask));
          outp.fc_sim_env_mean(iies,iiis,iG,igain) = mean(outp.fc_sim_env(mask));
          %                     outp.fc_sim_all(:,:,iies,iiis,iG) = outp.fc_sim_tmp;
          %         outp.fc_sim_env(iies,iiis,iG,igain) = mean(outp.fc_sim_env_tmp(mask));
          %
          outp.Ies(iies) = out.Ie;
          outp.Iis(iiis) = out.Ii;
          
          % KURAMOTO
          outp.kuramoto_mean (iies,iiis,iG,igain) = mean(out.KOPmean);
          outp.kuramoto_std (iies,iiis,iG,igain)  = mean(out.KOPsd);
          
          fclose all;
        catch me
          cnt = cnt + 1;
        end
        
        %
      end
    end
  end
end




save(sprintf('~/pmod/proc/pmod_wholebrain_woolrich_all_v%d.mat',v_sim),'outp')
%
else
    load(sprintf('~/pmod/proc/pmod_wholebrain_woolrich_all_v%d.mat',v_sim))
end

error('!')

%% PLOT BASIC PARAMETERS
% Lambda
% Peak frequency
% Correlation model with MEG: Lambda and FC
% ------------------------------------------
if (v==2 | v==3) && size(osc1,1)>26
  osc1 = osc1(11:36,1:26,:,:);
end

clear par
% ei = 2
idx = [find(round(Ies*100)/100==-3) find(round(Iis*100)/100==-3.5000)];
idx2 = [find(round(Ies*100)/100==-2) find(round(Iis*100)/100==-2.5000)];
Ies_old         = -4:0.1:-1;
Iis_old         = -5:0.1:-1;
igain = 1;
iG =1;
figure; set(gcf,'color','w')

% plot lambda
ax{1} = subplot(2,2,1); hold on
par = squeeze(mean(1./outp.lambda_env(:,:,:,iG)));
par(osc1(:,:,iG,igain)>0.5)=nan;
imagescnan(par,[100 200])
title('\lambda (model)');

% plot correlation FC model / MEG
ax{2} = subplot(2,2,2); hold on
par = squeeze(outp.r_env_rest_corr(:,:,iG));
par(osc1(:,:,iG,igain)>0.5)=nan;
imagescnan(par,[0 0.4])
title('Corr(FC_{sim}, FC_{exp})');

% plot correlation lambda model / MEG
ax{3} = subplot(2,2,3); hold on
par = squeeze(outp.lambda_r_rest(:,:,iG));
par(osc1(:,:,iG,igain)>0.5)=nan;
imagescnan(par,[0 0.15])
title('Corr(\lambda_{sim},\lambda_{exp})');

% plot peak freq model
ax{4} = subplot(2,2,4); hold on
par = squeeze(outp.peakfreq(:,:,iG));
par(osc1(:,:,iG,igain)>0.5)=nan;
imagescnan(par,[2 30])
title('Peak frequency');

for iax = 1 : length(ax)
  
  scatter(ax{iax},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
  scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
  
  line(ax{iax},[find(Iis==Iis_old(1)) find(Iis==Iis_old(1))],[find(Ies==Ies_old(1)) find(Ies==Ies_old(end))],'color','w')
  line(ax{iax},[find(Iis==Iis_old(1)) find(Iis==Iis_old(end))],[find(Ies==Ies_old(1)) find(Ies==Ies_old(1))],'color','w')
  line(ax{iax},[find(Iis==Iis_old(end)) find(Iis==Iis_old(end))],[find(Ies==Ies_old(1)) find(Ies==Ies_old(end))],'color','w')
  line(ax{iax},[find(Iis==Iis_old(end)) find(Iis==Iis_old(1))],[find(Ies==Ies_old(end)) find(Ies==Ies_old(end))],'color','w')
  
  
  
  ylabel(ax{iax},'Excitatory input'); xlabel(ax{iax},'Inhibitory input')
  set(ax{iax},'XTick',1:5:length(Iis),'XTickLabels',Iis(1:5:end))
  set(ax{iax},'YTick',1:5:length(Ies),'YTickLabels',Ies(1:5:end))
  
  colormap(plasma)
  tp_editplots(ax{iax})
end

% print(gcf,'-dpdf',sprintf('~/Dropbox/projects/phd/pmod/plots/pmod_wc_wholebrain_lambda_excitgain_iG%d_v%d.pdf',iG,v_sim))

%% PLOT BASIC PARAMETERS
% Lambda
% Peak frequency
% Correlation model with MEG: Lambda and FC
% ------------------------------------------
clear par
% ei = 2
idx = [find(round(Ies*100)/100==-3) find(round(Iis*100)/100==-3.5000)];
idx2 = [find(round(Ies*100)/100==-2) find(round(Iis*100)/100==-2.5000)];
Ies_old         = -4:0.1:-1;
Iis_old         = -5:0.1:-1;
igain = 1;
iG =1;
figure; set(gcf,'color','w')

% plot lambda
ax{1} = subplot(2,2,1); hold on
par = squeeze(outp.fc_sim_mean(:,:,iG,igain));
par(osc1(:,:,iG,igain)>0.5)=nan;
imagescnan(par,[0 0.05])
title('FCraw)');

% plot correlation FC model / MEG
ax{2} = subplot(2,2,2); hold on
par = squeeze(outp.fc_sim_env_mean(:,:,iG));
par(osc1(:,:,iG,igain)>0.5)=nan;
imagescnan(par,[0 0.05])
title('FCenv)');

% plot correlation FC model / MEG
ax{3} = subplot(2,2,3); hold on
par = squeeze(outp.r_rest_corr(:,:,iG,igain));
par(osc1(:,:,iG,igain)>0.5)=nan;
imagescnan(par,[0 0.3])
title('Corr(FCraw)');

% plot correlation FC model / MEG
ax{4} = subplot(2,2,4); hold on
par = squeeze(outp.r_env_rest_corr(:,:,iG,igain));
par(osc1(:,:,iG,igain)>0.5)=nan;
imagescnan(par,[0 0.3])
title('Corr(FCenv)');

for iax = 1 : length(ax)
  
  scatter(ax{iax},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
  scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
  
  line(ax{iax},[find(Iis==Iis_old(1)) find(Iis==Iis_old(1))],[find(Ies==Ies_old(1)) find(Ies==Ies_old(end))],'color','w')
  line(ax{iax},[find(Iis==Iis_old(1)) find(Iis==Iis_old(end))],[find(Ies==Ies_old(1)) find(Ies==Ies_old(1))],'color','w')
  line(ax{iax},[find(Iis==Iis_old(end)) find(Iis==Iis_old(end))],[find(Ies==Ies_old(1)) find(Ies==Ies_old(end))],'color','w')
  line(ax{iax},[find(Iis==Iis_old(end)) find(Iis==Iis_old(1))],[find(Ies==Ies_old(end)) find(Ies==Ies_old(end))],'color','w')
  
  
  
  ylabel(ax{iax},'Excitatory input'); xlabel(ax{iax},'Inhibitory input')
  set(ax{iax},'XTick',1:5:length(Iis),'XTickLabels',Iis(1:5:end))
  set(ax{iax},'YTick',1:5:length(Ies),'YTickLabels',Ies(1:5:end))
  
  colormap(plasma)
  tp_editplots(ax{iax})
end

% print(gcf,'-dpdf',sprintf('~/Dropbox/projects/phd/pmod/plots/pmod_wc_wholebrain_lambda_excitgain_iG%d_v%d.pdf',iG,v_sim))
