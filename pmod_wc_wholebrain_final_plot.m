%% FITTING
% pmod_wc_wholebrain_final_plot.m
%-------------------------------------------------------------------------
% VERSION 1: After meeting with tobi, 15-08-2018
%-------------------------------------------------------------------------
v           = 4;
Ies         = -4:0.1:-1;
Iis         = -5:0.1:-1;
Gg          = 0:0.2:1;
Gains       = 0:0.05:0.2;
nTrials     = 3;
tmax        = 6500; % in units of tauE
%-------------------------------------------------------------------------

v_sim = v;
% connectivity, AAL
v_conn =  1;
% simulations
% dfa, aal, lcmv
v_dfa = 2;

set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})
Gs = 1;

% ---------
% LOAD EMPIRICAL DFA (LCMV)
% ---------
load(sprintf(['~/pupmod/proc/conn/' 'pupmod_src_dfa_v%d.mat'],v_dfa));
% ---------
% MATCH DFA WITH AAL ORDER USED FOR SIMULATIONS
% ---------
pars = [];
pars.grid = 'medium';
pars.N = 90;
% ---------

dfa_emp_rest = nanmean(outp.dfa_all(:,:,1,1,1),2);
dfa_emp_task = nanmean(outp.dfa_all(:,:,1,1,2),2);

lambda_emp_rest = nanmean(outp.lambda_all(:,:,1,1,1),2);
lambda_emp_task = nanmean(outp.lambda_all(:,:,1,1,2),2);

load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v_conn));
load /home/tpfeffer/pupmod/proc/pow/pupmod_src_peakfreq_v3.mat
peakfreq_rest = m_res(1);
peakfreq_task = m_tsk(1);

load ~/pupmod/proc/conn/pupmod_all_kuramoto_v1.mat

mask = logical(tril(ones(90,90),-1));
mask = find(triu(ones(90))-eye(90));

fc_rest     =  squeeze(nanmean(cleandat(:,:,:,1,1,6),3));
fc_task     =  squeeze(nanmean(cleandat(:,:,:,1,2,6),3));
fc_rest_var =  std(nanmean(squeeze(nanmean(cleandat(:,:,:,1,1,6),3)))./max(nanmean(squeeze(nanmean(cleandat(:,:,:,1,1,6),3)))));

% if ~exist(sprintf('~/pmod/proc/pmod_wholebrain_rest_all_v%d.mat',v_sim))

for iies = 1 : length(Ies)
  iies
  for iiis = 1 : length(Iis)
    for iG = 1:length(Gg)
      for igain = 1:length(Gains)
        %         igain
        load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',iies,iiis,iG,igain,v_sim))

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
        
%         [~,peak_idx]=max(smooth(mean(out.PSD(out.f>3,:),2),20));
        outp.peakfreq(iies,iiis,iG,igain) = out.peakfreq;
        
        % get peak frequency from MEG data
%         outp.peakfreq_diff_res(iies,iiis,iG,igain) = outp.peakfreq(iies,iiis,iG,igain)-peakfreq_rest;
%         outp.peakfreq_diff_tsk(iies,iiis,iG,igain) = outp.peakfreq(iies,iiis,iG,igain)-peakfreq_task;
        
        pars.dim = 2;
        
        outp.fc_sim_tmp = tp_match_aal(pars,out.FC,3);  
        outp.fc_sim_env_tmp = tp_match_aal(pars,out.FC_env,3);  
        
%         outp.fc_sim_var(:,iies,iiis,iG,igain) = std(nanmean(outp.fc_sim_tmp)./max(nanmean(outp.fc_sim_tmp)));
        
        [outp.r_rest_corr(iies,iiis,iG,igain), outp.p_rest_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_tmp(mask),fc_rest(mask));
        [outp.r_env_rest_corr(iies,iiis,iG,igain), outp.p_env_rest_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_rest(mask));

%         [outp.r_rest_corr_avg(iies,iiis,iG,igain), outp.p_rest_corr_avg(iies,iiis,iG,igain)]=corr(nanmean(outp.fc_sim_tmp)',nanmean(fc_rest)');
        
%         outp.r_rest_corr_unc(iies,iiis,iG,igain) = dot(outp.fc_sim_tmp(mask),fc_rest(mask)) / sqrt(dot(outp.fc_sim_tmp(mask),outp.fc_sim_tmp(mask)) * dot(fc_rest(mask),fc_rest(mask)));
        %
        pars.dim = 1;
        
%         [outp.dfa_r_rest(iies,iiis,iG,igain), outp.dfa_p_rest(iies,iiis,iG,igain)] = corr(dfa_emp_rest(:),tp_match_aal(pars,repmat(outp.dfa_sim(:,iies,iiis,iG,igain),[1 90]),pars));
%        	[outp.dfa_env_r_rest(iies,iiis,iG,igain), outp.dfa_env_p_rest(iies,iiis,iG,igain)] = corr(dfa_emp_rest(:),tp_match_aal(pars,repmat(outp.dfa_env_sim(:,iies,iiis,iG,igain),[1 90]),pars));

        [outp.lambda_r_rest(iies,iiis,iG,igain), outp.lambda_p_rest(iies,iiis,iG,igain)] = corr(lambda_emp_rest,tp_match_aal(pars,repmat(outp.lambda(:,iies,iiis,iG,igain),[1 90]),pars));
        [outp.lambda_r_task(iies,iiis,iG,igain), outp.lambda_p_task(iies,iiis,iG,igain)] = corr(lambda_emp_task,tp_match_aal(pars,repmat(outp.lambda(:,iies,iiis,iG,igain),[1 90]),pars));
        
        [outp.lambda_env_r_rest(iies,iiis,iG,igain), outp.lambda_env_p_rest(iies,iiis,iG,igain)] = corr(lambda_emp_rest,tp_match_aal(pars,repmat(outp.lambda_env(:,iies,iiis,iG,igain),[1 90]),pars));

%         outp.dist_fc_rest (iies,iiis,iG,igain)  = mean(outp.fc_sim_tmp(mask))-mean(fc_rest(mask));
%         outp.dist_fc_task (iies,iiis,iG,igain)  = mean(outp.fc_sim_tmp(mask))-mean(fc_task(mask));
        
        outp.fc_sim_mean(iies,iiis,iG,igain) = mean(outp.fc_sim_tmp(mask));
%         outp.fc_sim_all(:,:,iies,iiis,iG,igain) = outp.fc_sim_tmp;
        outp.fc_sim_env_mean(iies,iiis,iG,igain) = mean(outp.fc_sim_env_tmp(mask));
        %
        outp.Ies(iies) = out.Ie;
        outp.Iis(iiis) = out.Ii;
        
        % KURAMOTO
        outp.kuramoto_mean (iies,iiis,iG,igain) = mean(out.KOPmean);
        outp.kuramoto_std (iies,iiis,iG,igain)  = mean(out.KOPsd);
        
%         outp.kuramoto_mean_diff (iies,iiis,iG,igain) = mean(out.KOPmean) - mean(kura_mean(:,1,1,1));
%         outp.kuramoto_std_diff (iies,iiis,iG,igain)  = mean(out.KOPsd)- mean(kura_std(:,1,1,1));
        
%         outp.psslp(:,iies,iiis,iG,igain)      = out.psslope;
%         outp.psslp_env(:,iies,iiis,iG,igain)  = out.psslope_env;
        
        fclose all;
        %
      end
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
% 
save(sprintf('~/pmod/proc/pmod_wholebrain_rest_all_v%d.mat',v_sim),'outp')
% % % 
% else
%   load(sprintf('~/pmod/proc/pmod_wholebrain_rest_all_v%d.mat',v_sim))
% end

error('!')

%% PLOT BASIC PARAMETERS: Functional connectivity
% ------------------------------------------
clear par

idx = [find(round(Ies*100)/100==-2.8) find(round(Iis*100)/100==-3.5000)];
idx2 = [find(round(Ies*100)/100==-1.8) find(round(Iis*100)/100==-2.4000)];

igain = 5;
G = 2;

figure; set(gcf,'color','w')

% plot FC
ax{1} = subplot(2,2,1); hold on
par = squeeze(outp.fc_sim_mean(:,:,G,igain));
% par(osc1>0.5)=nan;
imagescnan(par,[0 0.01])
title('FC_{FR}');

% plot FC env
ax{2} = subplot(2,2,2); hold on
par = squeeze(outp.fc_sim_env_mean(:,:,G,igain));
% par(osc1>0.5)=nan;
imagescnan(par,[0 0.01])
title('FC_{Env}');

% plot correlation FC sim w exp
ax{3} = subplot(2,2,3); hold on
par = squeeze(outp.r_rest_corr(:,:,G,igain));
% par(osc1>0.5)=nan;
imagescnan(par,[0 0.3])
title('r(FC_{FR})');

% plot peak freq model
ax{4} = subplot(2,2,4); hold on
par = squeeze(outp.r_env_rest_corr(:,:,G,igain));
% par(osc1>0.5)=nan;
imagescnan(par,[0 0.3])
title('r(FC_{Env})');

for iax = 1 : length(ax)
  scatter(ax{iax},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
  scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
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
  colormap(plasma)
  c = colorbar(ax{iax}); axis(ax{iax},'tight')
  c.Ticks = [min(c.Ticks) max(c.Ticks)];
 
end


print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_fc_gain%d_G%d_v%d.pdf',igain,G,v_sim))


%% PLOT BASIC PARAMETERS: Timescales
% ------------------------------------------
clear par

idx = [find(round(Ies*100)/100==-2.8) find(round(Iis*100)/100==-3.5000)];
idx2 = [find(round(Ies*100)/100==-1.8) find(round(Iis*100)/100==-2.4000)];

igain = 5
G = 1;

figure; set(gcf,'color','w')

% plot lambda
ax{1} = subplot(2,2,1); hold on
par = squeeze(mean(1./outp.lambda(:,:,:,G,igain)));
par(osc1>0.5)=nan;
imagescnan(par,[7 15])
title('\lambda_{FR}');

% plot correlation FC model / MEG
ax{2} = subplot(2,2,2); hold on
par = squeeze(mean(1./outp.lambda_env(:,:,:,G,igain)));
par(osc1>0.5)=nan;
imagescnan(par,[140 150])
title('\lambda_{Env}');

% plot correlation lambda model / MEG
ax{3} = subplot(2,2,3); hold on
par = squeeze(outp.lambda_r_rest(:,:,G,igain));
par(osc1>0.5)=nan;
imagescnan(par,[0 0.10])
title('r(\lambda_{FR})');

% plot peak freq model
ax{4} = subplot(2,2,4); hold on
par = squeeze(outp.lambda__env_r_rest(:,:,G,igain));
par(osc1>0.5)=nan;
imagescnan(par,[0 0.1])
title('r(\lambda_{Env})');


for iax = 1 : length(ax)
  scatter(ax{iax},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
  scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
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
  colormap(plasma)
  c = colorbar(ax{iax}); axis(ax{iax},'tight')
  c.Ticks = [min(c.Limits) max(c.Limits)];
 
end

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_lambda_gain%d_G%d_v%d.pdf',igain,G,v_sim))

%% PLOT BASIC PARAMETERS: DFA
% ------------------------------------------
clear par

idx = [find(round(Ies*100)/100==-2.8) find(round(Iis*100)/100==-3.5000)];
idx2 = [find(round(Ies*100)/100==-1.8) find(round(Iis*100)/100==-2.4000)];

igain = 1;
G = 1;

figure; set(gcf,'color','w')

% plot lambda
ax{1} = subplot(2,2,1); hold on
par = squeeze(mean(outp.dfa_sim(:,:,:,G,igain)));
par(osc1>0.5)=nan;
imagescnan(par,[0.5 0.7])
title('DFA_{FR}');

% plot correlation FC model / MEG
ax{2} = subplot(2,2,2); hold on
par = squeeze(mean(outp.dfa_env_sim(:,:,:,G,igain)));
par(osc1>0.5)=nan;
imagescnan(par,[0.50 0.70])
title('DFA_{Env}');

% plot correlation lambda model / MEG
ax{3} = subplot(2,2,3); hold on
par = squeeze(outp.dfa_r_rest(:,:,G,igain));
par(osc1>0.5)=nan;
imagescnan(par,[-0.2 0.20])
title('r(DFA_{FR})');

% plot peak freq model
ax{4} = subplot(2,2,4); hold on
par = squeeze(outp.r_env_rest_corr(:,:,G,igain));
par(osc1>0.5)=nan;
imagescnan(par,[0 0.2])
title('r(DFA_{Env})');


for iax = 1 : length(ax)
 scatter(ax{iax},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
  scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
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
  colormap(plasma)
  c = colorbar(ax{iax}); axis(ax{iax},'tight')
  c.Ticks = [min(c.Limits) max(c.Limits)];
 
end

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_dfa_gain%d_G%d_v%d.pdf',igain,G,v_sim))

%% PLOT GAIN VS NO GAIN
% increase in gain vs. baseline
clear par

idx = [find(round(Ies*100)/100==-2.8) find(round(Iis*100)/100==-3.5000)];
idx2 = [find(round(Ies*100)/100==-1.8) find(round(Iis*100)/100==-2.4000)];

igain = 5;
G = 1;

cmap = cbrewer('div', 'RdBu', 100,'pchip');
cmap = cmap(end:-1:1,:);

figure; set(gcf,'color','w')

% plot lambda
ax{1} = subplot(2,2,1); hold on
par = squeeze(abs(outp.fc_sim_mean(:,:,G,igain)))-abs(squeeze(outp.fc_sim_mean(:,:,G,3)));
par(osc1>0.5)=nan;
imagescnan(par,[-0.02 0.02])
title('Contrast: FC_{FR}');

% plot correlation lambda model / MEG
ax{2} = subplot(2,2,2); hold on
par = squeeze(abs(outp.fc_sim_env_mean(:,:,G,igain)))-abs(squeeze(outp.fc_sim_env_mean(:,:,G,3)));
par(osc1>0.5)=nan;
imagescnan(par,[-0.02 0.02])
title('Contrast: FC_{env}');

ax{3} = subplot(2,2,3); hold on
par = squeeze(mean(1./outp.lambda(:,:,:,G,igain)))-squeeze(1./mean(outp.lambda(:,:,:,G,3)));
par(osc1>0.5)=nan;
imagescnan(par,[-5 5])
title('Contrast: Lambda_{FR}');

% plot peak freq model
ax{4} = subplot(2,2,4); hold on
par = squeeze(1./mean(outp.lambda_env(:,:,:,G,igain)))-squeeze(1./mean(outp.lambda_env(:,:,:,G,3)));
par(osc1>0.5)=nan;
imagescnan(par,[-60 60])
title('Contrast: Lambda_{Env}');

for iax = 1 : length(ax)
 scatter(ax{iax},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
  scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
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
  colormap(cmap)
  c = colorbar(ax{iax}); axis(ax{iax},'tight')
  c.Ticks = [min(c.Limits) max(c.Limits)];
 
end


print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_gainvsbaseline_gain%d_G%d_v%d.pdf',igain,G,v_sim))


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

