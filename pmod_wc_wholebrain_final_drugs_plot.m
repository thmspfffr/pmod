%% FITTING
% pmod_wc_wholebrain_final_drugs_plot.m
% ---------------------------------
% Plot results from function above.
% Takes input from pmod_final_fitting, where parameters for rest and task
% are determined. From here, excitatory and inhibitory synaptic weights are
% changed parametrically.

%-------------------------------------------------------------------------
% VERSION 4: After meeting with tobi, 15-08-2018
%-------------------------------------------------------------------------
% v_sim       = 4;
% load(sprintf('~/pmod/proc/pmod_wc_wholebrain_rest_v%d.mat',v_sim))
% Ies         = [par.rest(1) par.task(1)];
% Iis         = [par.rest(2) par.task(2)];
% E_gain      = 0.75:0.025:1.25;
% I_gain      = 0.75:0.025:1.25;
% Gg          = 0.60;
% Gains       = 0;
% nTrials     = 3;
% tmax        = 6500; % in units of tauE
% wins = [2 20];
%-------------------------------------------------------------------------
% VERSION 5 After meeting with tobi, 15-08-2018
%-------------------------------------------------------------------------
v_sim           = 5;
load(sprintf('~/pmod/proc/pmod_wc_wholebrain_rest_v%d.mat',4))
Ies         = [par.rest(1) par.task(1)];
Iis         = [par.rest(2) par.task(2)];
E_gain      = 0.85:0.01:1.15;
I_gain      = 1;
Gg          = 0.60;
Gains       = -0.15:0.01:0.15;
nTrials     = 5;
tmax        = 6500; % in units of tauE
wins = [2 20];
%----------
% VERSION 6 After meeting with tobi, 15-08-2018
%-------------------------------------------------------------------------
v_sim       = 6;
load(sprintf('~/pmod/proc/pmod_wc_wholebrain_rest_v%d.mat',4))
Ies         = [par.rest(1) par.task(1)];
Iis         = [par.rest(2) par.task(2)];
E_gain      = 0.5:0.02:1.5;
I_gain      = 1;
Gg          = 0.60;
Gains       = -0.5:0.02:0.5;
nTrials     = 5;
tmax        = 6500; % in units of tauE
wins = [2 20];
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
%-------------------------------------------------------------------------
% ---------
% LOAD EMPIRICAL DFA (LCMV)
% ---------
load(sprintf(['~/pupmod/proc/conn/' 'pupmod_src_dfa_v%d.mat'],2));
% ---------
% MATCH DFA WITH AAL ORDER USED FOR SIMULATIONS
% ---------
pars = [];
pars.grid = 'medium';
pars.N = 90;
% ---------

dfa_emp_rest    = nanmean(outp.dfa_all(:,:,1,1,1),2);
lambda_emp_rest = nanmean(outp.lambda_all(:,:,1,1,1),2);

load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',1));
load /home/tpfeffer/pupmod/proc/pow/pupmod_src_peakfreq_v3.mat
load ~/pupmod/proc/conn/pupmod_all_kuramoto_v1.mat

peakfreq_rest = m_res(1);
peakfreq_task = m_tsk(1);
mask = find(triu(ones(90))-eye(90));

fc_rest     =  squeeze(nanmean(cleandat(:,:,:,1,1,6),3));

if ~exist(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_drugs_all_v%d.mat',v_sim))
  
  for exc_gain = 1 : length(E_gain)
    exc_gain
    for inh_gain = 1 : length(I_gain)
      inh_gain
      for iies = 1: length(Ies)
        for igain = 1 : length(Gains)
          
          % Load simulated output
          load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_drugs_Ie%d_ExcGain%d_InhGain%d_gain%d_v%d.mat',iies,exc_gain,inh_gain,igain,v_sim))
          % -----------
          % Time scales
          outp.lambda(:,exc_gain,inh_gain,igain,iies)      = mean(out.lambda,2);
          outp.lambda_env(:,exc_gain,inh_gain,igain,iies)  = mean(out.lambda_env,2);
          % DFA
          outp.dfa_sim(:,exc_gain,inh_gain,igain,iies)     = squeeze(mean(out.dfa,1));
          outp.dfa_env_sim(:,exc_gain,inh_gain,igain,iies) = squeeze(mean(out.dfa_env,1));
          % Peak freq
          outp.peakfreq(exc_gain,inh_gain,igain,iies) = out.peakfreq;
          % Simulated FC
          pars.dim = 2;
          outp.fc_sim_tmp = tp_match_aal(pars,out.FC,3);
          outp.fc_sim_env_tmp = tp_match_aal(pars,out.FC_env,3);
          outp.fc_sim_mean(exc_gain,inh_gain,igain,iies) = mean(outp.fc_sim_tmp(mask));
          outp.fc_sim_env_mean(exc_gain,inh_gain,igain,iies) = mean(outp.fc_sim_env_tmp(mask));
          % Correlations with MEG FC
          [outp.r_rest_corr(exc_gain,inh_gain,igain,iies), outp.p_rest_corr(exc_gain,inh_gain,igain,iies)]=corr(outp.fc_sim_tmp(mask),fc_rest(mask));
          [outp.r_env_rest_corr(exc_gain,inh_gain,igain,iies), outp.p_env_rest_corr(exc_gain,inh_gain,igain,iies)]=corr(outp.fc_sim_env_tmp(mask),fc_rest(mask));
          % Uncentered correlation
          outp.r_rest_corr_unc(exc_gain,inh_gain,igain,iies) = dot(outp.fc_sim_tmp(mask),fc_rest(mask)) / sqrt(dot(outp.fc_sim_tmp(mask),outp.fc_sim_tmp(mask)) * dot(fc_rest(mask),fc_rest(mask)));
          % Correlation with DFA
          pars.dim = 1;
          [outp.dfa_r_rest(exc_gain,inh_gain,igain,iies), outp.dfa_p_rest(exc_gain,inh_gain,igain,iies)] = corr(dfa_emp_rest(:),tp_match_aal(pars,repmat(outp.dfa_sim(:,exc_gain,inh_gain,igain,iies),[1 90]),pars));
          [outp.dfa_env_r_rest(exc_gain,inh_gain,igain,iies), outp.dfa_env_p_rest(exc_gain,inh_gain,igain,iies)] = corr(dfa_emp_rest(:),tp_match_aal(pars,repmat(outp.dfa_env_sim(:,exc_gain,inh_gain,igain,iies),[1 90]),pars));
          % Correlation lambda
          [outp.lambda_r_rest(exc_gain,inh_gain,igain,iies), outp.lambda_p_rest(exc_gain,inh_gain,igain,iies)] = corr(lambda_emp_rest,tp_match_aal(pars,repmat(outp.lambda(:,exc_gain,inh_gain,igain,iies),[1 90]),pars));
          [outp.lambda_env_r_rest(exc_gain,inh_gain,igain,iies), outp.lambda_env_p_rest(exc_gain,inh_gain,igain,iies)] = corr(lambda_emp_rest,tp_match_aal(pars,repmat(outp.lambda_env(:,exc_gain,inh_gain,igain,iies),[1 90]),pars));
          % Kuramoto
          outp.kuramoto_mean (exc_gain,inh_gain,igain,iies) = mean(out.KOPmean);
          outp.kuramoto_std (exc_gain,inh_gain,igain,iies)  = mean(out.KOPsd);
          
          fclose all;
          %
        end
      end
    end
  end
  save(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_drugs_all_v%d.mat',v_sim),'outp')
else
  load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_drugs_all_v%d.mat',v_sim))
end

error('!')


%% PLOT BASIC PARAMETERS: Functional connectivity
% ------------------------------------------
clear par

figure; set(gcf,'color','w')

% plot FC
ax{1} = subplot(2,2,1); hold on
par = squeeze(abs(outp.fc_sim_mean(:,1,:,1)))-abs(squeeze(outp.fc_sim_mean(16,1,16,1)));
% par(osc>0.5)=nan;
imagescnan(par,[-0.01 0.01])
title('FC_{FR}');

% plot FC
ax{2} = subplot(2,2,2); hold on
par = squeeze(abs(outp.fc_sim_mean(:,1,:,2)))-abs(squeeze(outp.fc_sim_mean(16,1,16,2)));
% par(osc>0.5)=nan;
imagescnan(par,[-0.01 0.01])
title('FC_{FR}');

% plot FC
ax{3} = subplot(2,2,3); hold on
par = squeeze(abs(outp.fc_sim_env_mean(:,1,:,1)))-abs(squeeze(outp.fc_sim_env_mean(16,1,16,1)));
% par(osc>0.5)=nan;
imagescnan(par,[-0.01 0.01])
title('FC_{Env}');

% plot FC
ax{4} = subplot(2,2,4); hold on
par = squeeze(abs(outp.fc_sim_env_mean(:,1,:,2)))-abs(squeeze(outp.fc_sim_env_mean(16,1,16,2)));
% par(osc>0.5)=nan;
imagescnan(par,[-0.01 0.01])
title('FC_{Env}');

cmap = cbrewer('div', 'RdBu', 100,'pchip'); cmap = cmap(end:-1:1,:);
for iax = 1 : length(ax)
%     scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
  if iax == 1
   	scatter(ax{iax},16,16,20,'markerfacecolor','w','markeredgecolor','k')
    ylabel(ax{iax},'E/I ratio');
    set(ax{iax},'YTick',1:10:length(I_gain ),'YTickLabels',num2cell(I_gain(1:10:end)))
    set(ax{iax},'XTick',1:10:length(E_gain),'XTickLabels',num2cell(E_gain(1:10:end)))
  elseif iax == 2
    scatter(ax{iax},16,16,20,'markerfacecolor','r','markeredgecolor','k')
    set(ax{iax},'YTick',1:10:length(I_gain ),'YTickLabels',num2cell(I_gain(1:10:end)))
    set(ax{iax},'XTick',1:10:length(E_gain),'XTickLabels',num2cell(E_gain(1:10:end)))
  elseif iax == 3
   	scatter(ax{iax},16,16,20,'markerfacecolor','w','markeredgecolor','k')
    xlabel(ax{iax},'Gain')
    ylabel(ax{iax},'E/I ratio');
    set(ax{iax},'YTick',1:10:length(I_gain ),'YTickLabels',num2cell(I_gain(1:10:end)))
    set(ax{iax},'XTick',1:10:length(E_gain),'XTickLabels',num2cell(E_gain(1:10:end)))
  elseif iax ==4
    scatter(ax{iax},16,16,20,'markerfacecolor','r','markeredgecolor','k')
    xlabel(ax{iax},'Gain')
    set(ax{iax},'YTick',1:10:length(I_gain ),'YTickLabels',num2cell(I_gain(1:10:end)))
    set(ax{iax},'XTick',1:10:length(E_gain),'XTickLabels',num2cell(E_gain(1:10:end)))
  end
  tp_editplots(ax{iax})
  colormap(ax{iax},redblue)
  c = colorbar(ax{iax}); axis(ax{iax},'tight')
  c.Ticks = c.Limits;
  
end


print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_final_drug_v%d.pdf',v_sim))


%% PLOT BASIC PARAMETERS: Timescales
% ------------------------------------------

figure; set(gcf,'color','w')

% plot FC
ax{1} = subplot(2,2,1); hold on
par = squeeze(abs(mean(1./outp.lambda(:,:,:,1))))-abs(squeeze(mean(1./outp.lambda(:,16,1,16,1))));
% par(osc>0.5)=nan;
imagescnan(par,[-7 7])
title('\lambda_{FR}');

% plot FC
ax{2} = subplot(2,2,2); hold on
par = squeeze(abs(mean(1./outp.lambda(:,:,:,2))))-abs(squeeze(mean(1./outp.lambda(:,11,11,:,2))));
% par(osc>0.5)=nan;
imagescnan(par,[-7 7])
title('\lambda_{FR}');

% plot FC
ax{3} = subplot(2,2,3); hold on
par = squeeze(abs(mean(1./outp.lambda_env(:,:,:,1))))-abs(squeeze(mean(1./outp.lambda_env(:,11,11,:,1))));
% par(osc>0.5)=nan;
imagescnan(par,[-7 7])
title('\lambda_{Env}');

% plot FC
ax{4} = subplot(2,2,4); hold on
par = squeeze(abs(mean(1./outp.lambda_env(:,:,:,2))))-abs(squeeze(mean(1./outp.lambda_env(:,11,11,:,2))));
% par(osc>0.5)=nan;
imagescnan(par,[-7 7])
title('\lambda_{Env}');

cmap = cbrewer('div', 'RdBu', 100,'pchip'); cmap = cmap(end:-1:1,:);
for iax = 1 : length(ax)
%     scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
  if iax == 1
   	scatter(ax{iax},11,11,20,'markerfacecolor','w','markeredgecolor','k')
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:10:length(I_gain ),'YTickLabels',num2cell(I_gain(1:10:end)))
    set(ax{iax},'XTick',1:10:length(E_gain),'XTickLabels',num2cell(E_gain(1:10:end)))
  elseif iax == 2
    scatter(ax{iax},11,11,20,'markerfacecolor','r','markeredgecolor','k')
    set(ax{iax},'YTick',1:10:length(I_gain ),'YTickLabels',num2cell(I_gain(1:10:end)))
    set(ax{iax},'XTick',1:10:length(E_gain),'XTickLabels',num2cell(E_gain(1:10:end)))
  elseif iax == 3
   	scatter(ax{iax},11,11,20,'markerfacecolor','w','markeredgecolor','k')
    xlabel(ax{iax},'Inhibitory input')
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:10:length(I_gain ),'YTickLabels',num2cell(I_gain(1:10:end)))
    set(ax{iax},'XTick',1:10:length(E_gain),'XTickLabels',num2cell(E_gain(1:10:end)))
  elseif iax ==4
    scatter(ax{iax},11,11,20,'markerfacecolor','r','markeredgecolor','k')
    xlabel(ax{iax},'Inhibitory input')
    set(ax{iax},'YTick',1:10:length(I_gain ),'YTickLabels',num2cell(I_gain(1:10:end)))
    set(ax{iax},'XTick',1:10:length(E_gain),'XTickLabels',num2cell(E_gain(1:10:end)))
  end
  tp_editplots(ax{iax})
  colormap(ax{iax},redblue)
  c = colorbar(ax{iax}); axis(ax{iax},'tight')
  c.Ticks = c.Limits;
  
end


print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_final_drug_lambda_v%d.pdf',v_sim))

%% PLOT BASIC PARAMETERS: DFA
% ------------------------------------------
figure; set(gcf,'color','w')

% plot FC
ax{1} = subplot(2,2,1); hold on
par = squeeze(abs(mean(outp.dfa_sim(:,:,:,1))))-abs(squeeze(mean(outp.dfa_sim(:,11,11,:,1))));
% par(osc>0.5)=nan;
imagescnan(par,[-0.05 0.05])
title('H_{FR}');

% plot FC
ax{2} = subplot(2,2,2); hold on
par = squeeze(abs(mean(outp.dfa_sim(:,:,:,2))))-abs(squeeze(mean(outp.dfa_sim(:,11,11,:,2))));
% par(osc>0.5)=nan;
imagescnan(par,[-0.05 0.05])
title('H_{FR}');

% plot FC
ax{3} = subplot(2,2,3); hold on
par = squeeze(abs(mean(outp.dfa_env_sim(:,:,:,1))))-abs(squeeze(mean(outp.dfa_env_sim(:,11,11,:,1))));
% par(osc>0.5)=nan;
imagescnan(par,[-0.05 0.05])
title('H_{Env}');

% plot FC
ax{4} = subplot(2,2,4); hold on
par = squeeze(abs(mean(outp.dfa_env_sim(:,:,:,2))))-abs(squeeze(mean(outp.dfa_env_sim(:,11,11,:,2))));
% par(osc>0.5)=nan;
imagescnan(par,[-0.05 0.05])
title('H_{Env}');

cmap = cbrewer('div', 'RdBu', 100,'pchip'); cmap = cmap(end:-1:1,:);
for iax = 1 : length(ax)
%     scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
  if iax == 1
   	scatter(ax{iax},11,11,20,'markerfacecolor','w','markeredgecolor','k')
    ylabel(ax{iax},'Inhibition gain');
    set(ax{iax},'YTick',1:10:length(I_gain ),'YTickLabels',num2cell(I_gain(1:10:end)))
    set(ax{iax},'XTick',1:10:length(E_gain),'XTickLabels',num2cell(E_gain(1:10:end)))
  elseif iax == 2
    scatter(ax{iax},11,11,20,'markerfacecolor','r','markeredgecolor','k')
    set(ax{iax},'YTick',1:10:length(I_gain ),'YTickLabels',num2cell(I_gain(1:10:end)))
    set(ax{iax},'XTick',1:10:length(E_gain),'XTickLabels',num2cell(E_gain(1:10:end)))
  elseif iax == 3
   	scatter(ax{iax},11,11,20,'markerfacecolor','w','markeredgecolor','k')
    xlabel(ax{iax},'Excitation gain')
    ylabel(ax{iax},'Inhibition gain');
    set(ax{iax},'YTick',1:10:length(I_gain ),'YTickLabels',num2cell(I_gain(1:10:end)))
    set(ax{iax},'XTick',1:10:length(E_gain),'XTickLabels',num2cell(E_gain(1:10:end)))
  elseif iax ==4
    scatter(ax{iax},11,11,20,'markerfacecolor','r','markeredgecolor','k')
    xlabel(ax{iax},'Excitation gain')
    set(ax{iax},'YTick',1:10:length(I_gain ),'YTickLabels',num2cell(I_gain(1:10:end)))
    set(ax{iax},'XTick',1:10:length(E_gain),'XTickLabels',num2cell(E_gain(1:10:end)))
  end
  tp_editplots(ax{iax})
  colormap(ax{iax},redblue)
  c = colorbar(ax{iax}); axis(ax{iax},'tight')
  c.Ticks = c.Limits;
  
end

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_final_drug_dfa_v%d.pdf',v_sim))


%% PLOT GAIN VS NO GAIN
% increase in gain vs. baseline
clear par

% idx = [find(round(Ies*100)/100==-2.8) find(round(Iis*100)/100==-3.5000)];
% idx2 = [find(round(Ies*100)/100==-1.8) find(round(Iis*100)/100==-2.4000)];

igain = 5;
G = 4;

cmap = cbrewer('div', 'RdBu', 256,'pchip');
cmap = cmap(end:-1:1,:);

figure; set(gcf,'color','w')

% plot lambda
ax{1} = subplot(2,2,1); hold on
par = squeeze(abs(outp.fc_sim_mean(:,:,G,igain)))-abs(squeeze(outp.fc_sim_mean(:,:,G,1)));
par(osc>0.5)=nan;
imagescnan(par,[-0.02 0.02])
title('Contrast: FC_{FR}');

% plot correlation lambda model / MEG
ax{2} = subplot(2,2,2); hold on
par = squeeze(abs(outp.fc_sim_env_mean(:,:,G,igain)))-abs(squeeze(outp.fc_sim_env_mean(:,:,G,1)));
par(osc>0.5)=nan;
imagescnan(par,[-0.05 0.05])
title('Contrast: FC_{env}');

ax{3} = subplot(2,2,3); hold on
par = squeeze(mean(1./outp.lambda(:,:,:,G,igain)))-squeeze(1./mean(outp.lambda(:,:,:,G,1)));
par(osc>0.5)=nan;
imagescnan(par,[-5 5])
title('Contrast: Lambda_{FR}');

% plot peak freq model
ax{4} = subplot(2,2,4); hold on
par = squeeze(1./mean(outp.lambda_env(:,:,:,G,igain)))-squeeze(1./mean(outp.lambda_env(:,:,:,G,1)));
par(osc>0.5)=nan;
imagescnan(par,[-30 30])
title('Contrast: Lambda_{Env}');

for iax = 1 : length(ax)
  %  scatter(ax{iax},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
  %   scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
  if iax == 1
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  elseif iax == 2
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  elseif iax == 3
    xlabel(ax{iax},'Inhibitory input')
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  elseif iax ==4
    xlabel(ax{iax},'Inhibitory input')
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  end
  tp_editplots(ax{iax})
  colormap(parula)
  c = colorbar(ax{iax}); axis(ax{iax},'tight')
  c.Ticks = [min(c.Limits) max(c.Limits)];
  
end


print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_final_gainvsbaseline_gain%d_G%d_v%d.pdf',igain,G,v_sim))


%%
clear par ax
figure; set(gcf,'color','w')
% clear

igain = 3;
% plot peak freq model
ax{1} = subplot(2,2,1); hold on
par = squeeze(outp.kuramoto_mean(:,:,G,igain));
par(osc1>0.5)=nan;
imagescnan(par,[0 0.3])
title('Kuramoto');

% plot peak freq model
ax{3} = subplot(2,2,3); hold on
par = squeeze(outp.kuramoto_mean(:,:,G,4))-squeeze(outp.kuramoto_mean(:,:,G,3));
par(osc1>0.5)=nan;
imagescnan(par,[-0.05 0.05])
title('Contrast: Kuramoto');

% plot peak freq model
ax{2} = subplot(2,2,2); hold on
par = squeeze(outp.kuramoto_std(:,:,G,igain));
par(osc1>0.5)=nan;
imagescnan(par,[0 0.0002])
title('Metastability');

% plot peak freq model
ax{4} = subplot(2,2,4); hold on
par = squeeze(outp.kuramoto_std(:,:,G,4))-squeeze(outp.kuramoto_std(:,:,G,3));
par(osc1>0.5)=nan;
imagescnan(par,[-0.0025 0.0025])
title('Contrast: Metastability');

for iax = 1 : length(ax)
  scatter(ax{iax},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
  scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
  if iax == 1
    colormap(ax{iax},plasma)
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  elseif iax == 2
    colormap(ax{iax},plasma)
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  elseif iax == 3
    colormap(ax{iax},cmap)
    xlabel(ax{iax},'Inhibitory input')
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  elseif iax ==4
    colormap(ax{iax},cmap)
    xlabel(ax{iax},'Inhibitory input')
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  end
  tp_editplots(ax{iax})
  
  c = colorbar(ax{iax}); axis(ax{iax},'tight')
  c.Ticks = [min(c.Limits) max(c.Limits)];
  
end


print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_kuramoto_gain%d_G%d_v%d.pdf',igain,G,v_sim))
%%
clear par ax
figure; set(gcf,'color','w')
% clear

igain = 3;
% plot peak freq model
ax{1} = subplot(2,2,1); hold on
par = squeeze(outp.degree(:,:,G,igain));
par(osc1>0.5)=nan;
imagescnan(par,[0.7 1.3])
title('Degree');

% plot peak freq model
ax{3} = subplot(2,2,3); hold on
par = squeeze(outp.degree(:,:,G,4))-squeeze(outp.degree(:,:,G,igain));
par(osc1>0.5)=nan;
imagescnan(par,[-0.3 0.3])
title('Contrast: Degree');

% plot peak freq
ax{2} = subplot(2,2,2); hold on
par = squeeze(outp.peakfreq(:,:,G,igain));
par(osc1>0.5)=nan;
imagescnan(par,[2 20])
title('Peak freq');

% plot peak freq model
ax{4} = subplot(2,2,4); hold on
par = squeeze(outp.peakfreq(:,:,G,4))-squeeze(outp.peakfreq(:,:,G,3));
par(osc1>0.5)=nan;
imagescnan(par,[-5 5])
title('Contrast: Freq');

for iax = 1 : length(ax)
  scatter(ax{iax},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
  scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
  if iax == 1
    colormap(ax{iax},plasma)
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  elseif iax == 2
    colormap(ax{iax},plasma)
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  elseif iax == 3
    colormap(ax{iax},cmap)
    xlabel(ax{iax},'Inhibitory input')
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  elseif iax ==4
    colormap(ax{iax},cmap)
    xlabel(ax{iax},'Inhibitory input')
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  end
  tp_editplots(ax{iax})
  
  c = colorbar(ax{iax}); axis(ax{iax},'tight')
  c.Ticks = [min(c.Limits) max(c.Limits)];
  
end

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_degreefreq_gain%d_G%d_v%d.pdf',igain,G,v_sim))
%%
cmap = cbrewer('div', 'RdBu', 100,'pchip');
cmap = cmap(end:-1:1,:);

clear par ax
figure; set(gcf,'color','w')
% clear

igain = 3;
% plot peak freq model
ax{1} = subplot(2,2,1); hold on
par = squeeze(outp.peakfreq_diff_res(:,:,1,igain));
par(osc1>0.5)=nan;
imagescnan(par,[-3 3])
title('Peak freq: Difference');


% plot peak freq model
ax{2} = subplot(2,2,2); hold on
par = squeeze(outp.peakfreq_diff_res(:,:,1,igain));
par(abs(par)<1)=1; par(abs(par)>1)=0;
par(osc1>0.5)=nan; m1 = par>0;
imagescnan(par,[-1 1])
title('Peak freq: Masked');

% plot peak freq model
ax{3} = subplot(2,2,3); hold on
par = squeeze(outp.r_env_rest_corr(:,:,1,igain));
par(osc1>0.5)=nan;
imagescnan(par,[0 0.3])
title('Peak freq: Masked');

igain = 3;
% plot peak freq model
ax{4} = subplot(2,2,4); hold on
par = squeeze(outp.r_env_rest_corr(:,:,1,igain));
par(abs(par)<0.2)=0; par(abs(par)>0.2)=1;
par(osc1>0.5)=nan; m2 = par>0;
imagescnan(par,[-1 1])
title('Peak freq: Masked');

for iax = 1 : 4
  scatter(ax{iax},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
  %   scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
  if iax == 1
    colormap(ax{iax},cmap)
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  elseif iax == 2
    colormap(ax{iax},cmap)
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  elseif iax == 3
    colormap(ax{iax},plasma)
    xlabel(ax{iax},'Inhibitory input')
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  elseif iax ==4
    colormap(ax{iax},cmap)
    xlabel(ax{iax},'Inhibitory input')
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  end
  tp_editplots(ax{iax})
  
  c = colorbar(ax{iax}); axis(ax{iax},'tight')
  c.Ticks = [min(c.Limits) max(c.Limits)];
  
end

figure; set(gcf,'color','w')
ax{1} = subplot(2,2,1); hold on
par = double(m1&m2);
par(osc1>0.5)=nan;
imagescnan(par,[-1 1])
title('Peak freq: Masked');
colormap(gca,cmap)

[i,k]=find(m1&m2)
idx = [round(mean(i)) round(mean(k))];

d_frest_ftask = peakfreq_task-peakfreq_rest;

% plot peak freq model
ax{2} = subplot(2,2,2); hold on
par = (squeeze(outp.peakfreq(:,:,1,igain))-outp.peakfreq(idx(1),idx(2),1,igain))-d_frest_ftask;
par(abs(par)<1)=1; par(abs(par)>1)=0;
par(osc1>0.5)=nan; m1 = par>0;
imagescnan(par,[-1 1])
title('Peak freq: Difference');

% plot peak freq model
ax{3} = subplot(2,2,3); hold on
par = double(m1&m2);
par(osc1>0.5)=nan;

imagescnan(par,[-1 1])
title('Peak freq: Difference');
colormap(cmap)


[i,k]=find(m1&m2)
idx2 = [round(mean(i)) round(mean(k))];

ax{4} = subplot(2,2,4); hold on
par =zeros(size((par)));
par(osc1>0.5)=nan;
imagescnan(par,[-1 1]);
title('Peak freq: Difference');


for iax = 1 : 4
  scatter(ax{iax},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
  if iax == 1
    colormap(ax{iax},plasma)
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  elseif iax == 2
    colormap(ax{iax},plasma)
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  elseif iax == 3
    colormap(ax{iax},plasma)
    xlabel(ax{iax},'Inhibitory input')
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  elseif iax == 4
    colormap(ax{iax},plasma)
    xlabel(ax{iax},'Inhibitory input')
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
    scatter(ax{4},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
    
  end
  tp_editplots(ax{iax})
  
  %   c = colorbar(ax{iax});
  axis(ax{iax},'tight')
  %   c.Ticks = [min(c.Limits) max(c.Limits)];
  %
end

clear par
par.rest = [Ies(idx(2)) Iis(idx(1))];
par.task = [Ies(idx2(2)) Iis(idx2(1))];
par.descr = 'First entry: Excitation (Ies), Second entry: Inhibition (Iis)';
save(sprintf('~/pmod/proc/pmod_wc_wholebrain_rest_v%d.mat',v_sim),'par')



%%
igain = 3;
% plot peak freq model
ax{1} = subplot(2,2,1); hold on
par = squeeze(outp.peakfreq_diff_res(:,:,1,igain));
par(osc1>0.5)=nan;
imagescnan(par,[-15 15])
title('Degree');
colormap(plasma)

scatter(gca,idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
scatter(gca,idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')

%% SPATIAL MAPS
addpath /home/gnolte/meg_toolbox/meg/

load /home/gnolte/meth/templates/sa_template
load /home/gnolte/meth/templates/mri.mat
sa_template.grid_cortex400 = select_chans(sa_template.grid_cortex3000,400);

grid  = sa_template.grid_cortex400;
g1    = sa_template.grid_cortex400;
g2    = sa_template.cortex10K.vc;
vc    = sa_template.vc;
dd    = .75;

d = mean(outp.fc_sim_all(:,:,idx(1),idx(2),1,3));

par_interp = spatfiltergauss(d,g1,dd,g2);

para =[];
para.colorlimits = [0 0.03];


% PLOT RESULTS
para.filename = sprintf('~/pmod/plots/all_src_tsk_v%d.png',v_sim);
tp_showsource(par_interp,cmap,sa_template,para);

