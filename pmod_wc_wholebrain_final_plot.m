%% FITTING
% pmod_wc_wholebrain_final_plot.m
%-------------------------------------------------------------------------
% VERSION 2: After meeting with tobi, 24-08-2018: even more fine gained
% %-------------------------------------------------------------------------
v           = 2;
Ies         = -4:0.025:-1;
Iis         = -5:0.025:-2;
Gg          = 0.6;
Gains       = 0;
nTrials     = 1;
tmax        = 6500; % in units of tauE
wins = [2 20]; 
%-------------------------------------------------------------------------
% VERSION 1: After meeting with tobi, 15-08-2018
%-------------------------------------------------------------------------
% v           = 4;
% Ies         = -4:0.1:-1;
% Iis         = -5:0.1:-1;
% Gg          = 0:0.2:1;
% Gains       = 0:0.05:0.2;
% nTrials     = 3;
% tmax        = 6500; % in units of tauE
%-------------------------------------------------------------------------
% v           = 3;
% Ies         = -4:0.01:-1;
% Iis         = -5:0.01:-1;
% Gg          = 0.6;
% Gains       = 0;
% nTrials     = 1;
% tmax        = 6500; % in units of tauE
% wins = [2 20]; 
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

clear outp

load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v_conn));
load /home/tpfeffer/pupmod/proc/pow/pupmod_src_peakfreq_v3.mat
peakfreq_rest = m_res(1);
peakfreq_task = m_tsk(1);

load ~/pupmod/proc/conn/pupmod_all_kuramoto_v1.mat

mask = logical(tril(ones(90,90),-1));

pars.dim = 2;
pars.grid = 'medium';
pars.N = 90;
pars.transfer = 'to_bcn';

fc_rest     =  tp_match_aal(pars,squeeze(nanmean(cleandat(:,:,:,1,1,6),3)));
fc_task     =  tp_match_aal(pars,squeeze(nanmean(cleandat(:,:,:,1,2,6),3)));
fc_rest_var =  std(nanmean(squeeze(nanmean(cleandat(:,:,:,1,1,6),3)))./max(nanmean(squeeze(nanmean(cleandat(:,:,:,1,1,6),3)))));

fc_rest_indiv = reshape(squeeze(cleandat(:,:,:,1,1,6)),[90*90 28]);
fc_rest_indiv = fc_rest_indiv(mask,:);
fc_task_indiv = reshape(squeeze(cleandat(:,:,:,1,2,6)),[90*90 28]);
fc_task_indiv = fc_task_indiv(mask,:);

if ~exist(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_all_v%d.mat',v_sim))

for iies = 1 : length(Ies)
  iies
  for iiis = 1 : length(Iis)
%     iiis
    for iG = 1:length(Gg)
      for igain = 1:length(Gains)
        %         igain
        load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',iies,iiis,iG,igain,v_sim))

        % Time scales
%         outp.lambda(:,iies,iiis,iG,igain)      = mean(out.lambda,2);
%         outp.lambda_env(:,iies,iiis,iG,igain)  = mean(out.lambda_env,2);
        % DFA
%         outp.dfa_sim(:,iies,iiis,iG,igain)     = squeeze(mean(out.dfa,1));
%         outp.dfa_env_sim(:,iies,iiis,iG,igain) = squeeze(mean(out.dfa_env,1));
        
%         [~,peak_idx]=max(smooth(mean(out.PSD(out.f>3,:),2),20));
        outp.peakfreq(iies,iiis,iG,igain) = out.peakfreq;
        
        pars.dim = 2;
        
        outp.fc_sim_tmp = out.FC;  
        outp.fc_sim_env_tmp = out.FC_env;  
        
%         outp.fc_sim_tmp = tp_match_aal(pars,out.FC,3);  
%         outp.fc_sim_env_tmp = tp_match_aal(pars,out.FC_env,3);  
                
        [outp.r_rest_corr(iies,iiis,iG,igain), outp.p_rest_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_tmp(mask),fc_rest(mask));
        [outp.r_env_rest_corr(iies,iiis,iG,igain), outp.p_env_rest_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_rest(mask));
        [outp.r_env_task_corr(iies,iiis,iG,igain), outp.p_env_task_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_task(mask));
        [outp.r_env_rest_indiv_corr(:,iies,iiis,iG,igain), outp.p_env_rest_indiv_corr(:,iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_rest_indiv);

%         [outp.r_rest_corr_avg(iies,iiis,iG,igain), outp.p_rest_corr_avg(iies,iiis,iG,igain)]=corr(nanmean(outp.fc_sim_tmp)',nanmean(fc_rest)');
%         outp.r_rest_corr_unc(iies,iiis,iG,igain) = dot(outp.fc_sim_tmp(mask),fc_rest(mask)) / sqrt(dot(outp.fc_sim_tmp(mask),outp.fc_sim_tmp(mask)) * dot(fc_rest(mask),fc_rest(mask)));
        %
        pars.dim = 1;
        
%         [outp.dfa_r_rest(iies,iiis,iG,igain), outp.dfa_p_rest(iies,iiis,iG,igain)] = corr(dfa_emp_rest(:),tp_match_aal(pars,repmat(outp.dfa_sim(:,iies,iiis,iG,igain),[1 90]),pars));
%        	[outp.dfa_env_r_rest(iies,iiis,iG,igain), outp.dfa_env_p_rest(iies,iiis,iG,igain)] = corr(dfa_emp_rest(:),tp_match_aal(pars,repmat(outp.dfa_env_sim(:,iies,iiis,iG,igain),[1 90]),pars));

%         [outp.lambda_r_rest(iies,iiis,iG,igain), outp.lambda_p_rest(iies,iiis,iG,igain)] = corr(lambda_emp_rest,tp_match_aal(pars,repmat(outp.lambda(:,iies,iiis,iG,igain),[1 90]),pars));
%         [outp.lambda_r_task(iies,iiis,iG,igain), outp.lambda_p_task(iies,iiis,iG,igain)] = corr(lambda_emp_task,tp_match_aal(pars,repmat(outp.lambda(:,iies,iiis,iG,igain),[1 90]),pars));
        
%         [outp.lambda_env_r_rest(iies,iiis,iG,igain), outp.lambda_env_p_rest(iies,iiis,iG,igain)] = corr(lambda_emp_rest,tp_match_aal(pars,repmat(outp.lambda_env(:,iies,iiis,iG,igain),[1 90]),pars));
        
        outp.fc_sim_mean(iies,iiis,iG,igain) = mean(outp.fc_sim_tmp(mask));
        outp.fc_sim_env_mean(iies,iiis,iG,igain) = mean(outp.fc_sim_env_tmp(mask));
        %
        outp.Ies(iies) = out.Ie;
        outp.Iis(iiis) = out.Ii;
        
        % KURAMOTO
        outp.kuramoto_mean (iies,iiis,iG,igain) = mean(out.KOPmean);
        outp.kuramoto_std (iies,iiis,iG,igain)  = mean(out.KOPsd);
                
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
  save(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_all_v%d.mat',v_sim),'outp')
% % 
else
  load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_all_v%d.mat',v_sim))
end

error('!')


%% PLOT BASIC PARAMETERS: Functional connectivity
% ------------------------------------------
clear par
% osc = osc1(:,:,4);
nancol = [0.97 0.97 0.97];
% idx = [find(round(Ies*100)/100==-2.8) find(round(Iis*100)/100==-3.5000)];
% idx2 = [find(round(Ies*100)/100==-1.8) find(round(Iis*100)/100==-2.4000)];

igain = 1;
G = 4;
oscthres = 0;
figure; set(gcf,'color','w')

% plot FC
ax{1} = subplot(2,2,1); hold on
par = squeeze(outp.fc_sim_mean(:,:,G,igain));
par(osc>oscthres)=nan;
imagescnan(par,[0 0.02],'NanColor',nancol)
title('FC_{FR}');

% plot FC env
ax{2} = subplot(2,2,2); hold on
par = squeeze(outp.fc_sim_env_mean(:,:,G,igain));
par(osc>oscthres)=nan;
imagescnan(par,[0 0.02],'NanColor',nancol)
title('FC_{Env}');

% plot correlation FC sim w exp
ax{3} = subplot(2,2,3); hold on
par = squeeze(outp.r_rest_corr(:,:,G,igain));
par(osc>oscthres)=nan;
imagescnan(par,[0 0.3],'NanColor',nancol)
title('r(FC_{FR})');

% plot peak freq model
ax{4} = subplot(2,2,4); hold on
par = squeeze(outp.r_env_rest_corr(:,:,G,igain));
par(osc>oscthres)=nan;
imagescnan(par,[0 0.3],'NanColor',nancol)
title('r(FC_{Env})');

for iax = 1 : length(ax)
%   scatter(ax{4},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
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
  colormap(plasma)
  c = colorbar(ax{iax}); axis(ax{iax},'equal');
  axis(ax{iax},[1 length(Iis) 1 length(Ies) ])
  c.Ticks = [min(c.Ticks) max(c.Ticks)];
   axis(ax{iax},'square')
 
end

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_final_fc_gain%d_G%d_v%d.pdf',igain,G,v_sim))

%% PLOT BASIC PARAMETERS: Timescales
% ------------------------------------------
clear par
oscthresh = 0;
% idx = [find(round(Ies*100)/100==-2.8) find(round(Iis*100)/100==-3.5000)];
% idx2 = [find(round(Ies*100)/100==-1.8) find(round(Iis*100)/100==-2.4000)];

igain = 1
G = 4;

figure; set(gcf,'color','w')

% plot lambda
ax{1} = subplot(2,2,1); hold on
par = squeeze(mean(1./outp.lambda(:,:,:,G,igain)));
par(osc>oscthresh)=nan;
imagescnan(par,[7 15])
title('\lambda_{FR}');

% plot correlation FC model / MEG
ax{2} = subplot(2,2,2); hold on
par = squeeze(mean(1./outp.lambda_env(:,:,:,G,igain)));
par(osc>oscthresh)=nan;
imagescnan(par,[135 145])
title('\lambda_{Env}');

% plot correlation lambda model / MEG
ax{3} = subplot(2,2,3); hold on
par = squeeze(outp.lambda_r_rest(:,:,G,igain));
par(osc>oscthresh)=nan;
imagescnan(par,[0 0.10])
title('r(\lambda_{FR})');

% plot peak freq model
ax{4} = subplot(2,2,4); hold on
par = squeeze(outp.lambda_env_r_rest(:,:,G,igain));
par(osc>oscthresh)=nan;
imagescnan(par,[0 0.1])
title('r(\lambda_{Env})');


for iax = 1 : length(ax)
%   scatter(ax{iax},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
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
  colormap(plasma)
  c = colorbar(ax{iax}); axis(ax{iax},'tight')
  c.Ticks = [min(c.Limits) max(c.Limits)];
 
end

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_final_lambda_gain%d_G%d_v%d.pdf',igain,G,v_sim))

%% PLOT GAIN VS NO GAIN
% increase in gain vs. baseline
clear par
load(sprintf('~/pmod/proc/pmod_final_fitting_indivfits_taskandrest_v%d.mat',v_sim))

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
par(osc>oscthres)=nan;
imagescnan(par,[-0.02 0.02],'NanColor',nancol)
title('Contrast: FC_{FR}');

% plot correlation lambda model / MEG
ax{2} = subplot(2,2,2); hold on
par = squeeze(abs(outp.fc_sim_env_mean(:,:,G,igain)))-abs(squeeze(outp.fc_sim_env_mean(:,:,G,1)));
par(osc>oscthres)=nan;
imagescnan(par,[-0.02 0.02],'NanColor',nancol)
title('Contrast: FC_{env}');

ax{3} = subplot(2,2,3); hold on
par = squeeze(mean(1./outp.lambda(:,:,:,G,igain)))-squeeze(1./mean(outp.lambda(:,:,:,G,1)));
par(osc>oscthres)=nan;
imagescnan(par,[-5 5],'NanColor',nancol)
title('Contrast: Lambda_{FR}');

% plot peak freq model
ax{4} = subplot(2,2,4); hold on
par = squeeze(1./mean(outp.lambda_env(:,:,:,G,igain)))-squeeze(1./mean(outp.lambda_env(:,:,:,G,1)));
par(osc>oscthres)=nan;
imagescnan(par,[-30 30],'NanColor',nancol)
title('Contrast: Lambda_{Env}');

for iax = 1 : length(ax)
 scatter(ax{iax},indiv_idx.rest(:,2),indiv_idx.rest(:,1),20,'markerfacecolor','w','markeredgecolor','k')
  scatter(ax{iax},indiv_idx.task(:,2),indiv_idx.task(:,1),20,'markerfacecolor','y','markeredgecolor','k')
  c = colorbar(ax{iax}); 
 c.Ticks = [c.Limits];
%   scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
  if iax == 1
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:10:length(Ies ),'YTickLabels',num2cell(Ies(1:10:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
%     c.Position = [0.44 0.588 0.01 0.3310];
  elseif iax == 2
    set(ax{iax},'YTick',1:10:length(Ies ),'YTickLabels',num2cell(Ies(1:10:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
%     c.Position = [0.88 0.588 0.01 0.3310];
  elseif iax == 3
    xlabel(ax{iax},'Inhibitory input')
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:10:length(Ies ),'YTickLabels',num2cell(Ies(1:10:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
%     c.Position = [0.44 0.114 0.01 0.3310];
  elseif iax ==4
    xlabel(ax{iax},'Inhibitory input')
    set(ax{iax},'YTick',1:10:length(Ies ),'YTickLabels',num2cell(Ies(1:10:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
%     c.Position = [0.88 0.114 0.01 0.3310];
  end
  tp_editplots(ax{iax})
  colormap(redblue)
  axis(ax{iax},'tight')
  axis(ax{iax},[1 length(Iis) 1 length(Ies) ])
 
   axis(ax{iax},'square')
   axis(ax{iax},'on'); 
end

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_final_gainvsbaseline_gain%d_G%d_v%d.pdf',igain,G,v_sim))

indiv_idx.task(indiv_idx.rest(:,1)==0,:) = [];
indiv_idx.rest(indiv_idx.rest(:,1)==0,:) = [];

tmp1 = squeeze(abs(outp.fc_sim_env_mean(indiv_idx.rest(:,1),indiv_idx.rest(:,2),G,igain)))-abs(squeeze(outp.fc_sim_env_mean(indiv_idx.rest(:,1),indiv_idx.rest(:,2),:,G,1)));
diag_mask = eye(size(tmp1,1),size(tmp1,2)); 
rest = tmp1(logical(diag_mask));

tmp1 = squeeze(abs(outp.fc_sim_env_mean(indiv_idx.task(:,1),indiv_idx.task(:,2),G,igain)))-abs(squeeze(outp.fc_sim_env_mean(indiv_idx.task(:,1),indiv_idx.task(:,2),:,G,1)));
diag_mask = eye(size(tmp1,1),size(tmp1,2)); 
task = tmp1(logical(diag_mask));

figure;set(gcf,'color','w')
subplot(2,2,1); hold on
scatter(ones(size(rest,1),1)-(0.5-rand(size(rest,1),1))/2,rest,30,'markeredgecolor','k','markerfacecolor','w')
scatter(2*ones(size(task,1),1)-(0.5-rand(size(rest,1),1))/2,task,30,'markeredgecolor','k','markerfacecolor','y')

axis([0 3 -0.2 1]); axis square; colorbar; tp_editplots
print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_final_gainvsbaseline_gain%d_bar_G%d_v%d.pdf',igain,G,v_sim))

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
igain = 1;
% plot peak freq model
ax{1} = subplot(2,2,1); hold on
par = squeeze(outp.fc_sim_env_mean(:,:,4,igain));
par(osc>0)=nan;
imagescnan(par,[0 0.02])
title('Degree');
colormap(plasma)

clear par
load(sprintf('~/pmod/proc/pmod_final_fitting_fits_v%d.mat',v_sim))

idx = [find(Ies==par.rest(1)) find(Iis==par.rest(2))];

scatter(gca,idx(1),idx(2),20,'markerfacecolor','w','markeredgecolor','k')

%% PARAMETERS OF INTEREST

load(sprintf('~/pmod/proc/pmod_final_fitting_indivfits_taskandrest_v%d.mat',v_sim))

% outp.fc_sim_mean(find(Ies==par.rest(1)),find(Iis==par.rest(2)),4,1)
% outp.fc_sim_mean(find(Ies==par.task(1)),find(Iis==par.task(2)),4,1)
% 
% figure; set(gcf,'color','w');
% 
% mean(outp.lambda_env(:,find(Ies==par.rest(1)),find(Iis==par.rest(2)),4,1))
% mean(outp.lambda_env(:,find(Ies==par.task(1)),find(Iis==par.task(2)),4,1))
% 
% 
% outp.peakfreq(find(Ies==par.rest(1)),find(Iis==par.rest(2)),4,1)
% outp.peakfreq(find(Ies==par.task(1)),find(Iis==par.task(2)),4,1)
% 
% % BAR PLOTS
diag_mask = logical(eye(22,22));
igain = 4

% get rid of subjects for which fit didn't work
nanidx = find(~isnan(indiv_idx.rest(:,1)));
indiv_idx.task=indiv_idx.task(nanidx,:);
indiv_idx.rest=indiv_idx.rest(nanidx,:);

par_rest = squeeze(abs(outp.fc_sim_env_mean(indiv_idx.rest(:,1),indiv_idx.rest(:,2),4,igain)))-abs(squeeze(outp.fc_sim_env_mean(indiv_idx.rest(:,1),indiv_idx.rest(:,2),4,1)))
par_rest = par_rest(diag_mask);

par_task = squeeze(abs(outp.fc_sim_env_mean(indiv_idx.task(:,1),indiv_idx.task(:,2),4,igain)))-abs(squeeze(outp.fc_sim_env_mean(indiv_idx.task(:,1),indiv_idx.task(:,2),4,1)))
par_task = par_task(diag_mask);


figure;set(gcf,'color','w')
subplot(2,2,1)
bar([1:22],[par_rest])
tp_editplots; axis square
axis([0 23 -0.2 1])
subplot(2,2,2)
bar([1:22],[par_task])
tp_editplots; axis square
axis([0 23 -0.1 1])

% print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_final_gainvsbaseline_bar_igain%d_v%d.pdf',igain,v_sim))

fc = cleandat(:,:,nanidx,2,2,6)-cleandat(:,:,nanidx,1,2,6);
fc = fc(repmat(mask,[1 1 22]));

%% FC RESULTS

M_rest = reshape(abs(cleandat(:,:,:,:,:,6)),[90*90 28 3 2]);
M_rest = M_rest(find(mask),:,:,:,:);
M_rest = squeeze(mean(M_rest,1));

delta_prc_all = 100* (M_rest(:,2,1)-M_rest(:,1,1)) ./ M_rest(:,1,1);
delta_prc_mean(1) = mean(delta_prc_all);
delta_prc_sem(1)  = std(delta_prc_all)/sqrt(28);

delta_prc_all = 100* (M_rest(:,2,2)-M_rest(:,1,2)) ./ M_rest(:,1,2);
delta_prc_mean(2) = mean(delta_prc_all);
delta_prc_sem(2)  = std(delta_prc_all)/sqrt(28);

figure; set(gcf,'color','white')
bar([1 2],delta_prc_mean);
axis([0 3 -2 18])

M_rest = reshape(abs(cleandat(:,:,:,:,:,6)),[90*90 28 3 2]);
M_rest = M_rest(find(mask),:,:,:,:);
M_rest = squeeze(mean(M_rest,1));

delta_prc_all = 100* (M_rest(:,3,1)-M_rest(:,1,1)) ./ M_rest(:,1,1);
delta_prc_mean(1) = mean(delta_prc_all);
delta_prc_sem(1)  = std(delta_prc_all)/sqrt(28);

delta_prc_all = 100* (M_rest(:,3,2)-M_rest(:,1,2)) ./ M_rest(:,1,2);
delta_prc_mean(2) = mean(delta_prc_all);
delta_prc_sem(2)  = std(delta_prc_all)/sqrt(28);

figure; set(gcf,'color','white')
bar([1 2],delta_prc_mean);
axis([0 3 -18 2])


%%

M_rest = reshape(abs(cleandat(:,:,:,:,:,6)),[90*90 28 3 2]);
M_rest = M_rest(find(mask),:,:,:,:);
M_rest = squeeze(mean(M_rest,1));

delta_prc_all = 100* (M_rest(:,2,1)-M_rest(:,1,1)) ./ M_rest(:,1,1);
delta_prc_mean(1) = mean(delta_prc_all);
delta_prc_sem(1)  = std(delta_prc_all)/sqrt(28);

delta_prc_all = 100* (M_rest(:,2,2)-M_rest(:,1,2)) ./ M_rest(:,1,2);
delta_prc_mean(2) = mean(delta_prc_all);
delta_prc_sem(2)  = std(delta_prc_all)/sqrt(28);

figure; set(gcf,'color','white')
bar([1 2],delta_prc_mean);
axis([0 3 -2 18])


M_rest = reshape(abs(cleandat(:,:,:,:,:,:)),[90*90 28 3 2 13]);
M_rest = M_rest(find(mask),:,:,:,:,:);
M_rest = squeeze(mean(M_rest,1));

[~,p]=ttest(M_rest(:,3,1,:),M_rest(:,1,1,:))
foi_range 	= unique(round(2.^[1:0.5:7]));

figure; set(gcf,'color','white'); hold on

ax{1} = subplot(2,2,1); hold on
plot(1:13,squeeze(p),'linewidth',3)
line([0 13],[0.05 0.05],'linestyle',':','color',[0.7 0.7 0.7])
plot([6 7],squeeze(p(6:7)),'k.','markersize',25)
axis([0 14 -0.1 1])

[~,p]=ttest(M_rest(:,3,2,:),M_rest(:,1,2,:))
plot(1:13,squeeze(p),'linewidth',3,'color',[0.4 0.8 0.9])
line([0 13],[0.05 0.05],'linestyle',':','color',[0.7 0.7 0.7])
% plot([6 7],squeeze(p(6:7)),'k.','markersize',25)
axis([0 14 -0.1 1])

ax{2} = subplot(2,2,2); hold on
[~,p]=ttest(M_rest(:,3,2,:)-M_rest(:,1,2,:),M_rest(:,3,1,:)-M_rest(:,1,1,:))
plot(1:13,squeeze(p),'linewidth',3,'color','k')
line([0 13],[0.05 0.05],'linestyle',':','color',[0.7 0.7 0.7])
% plot([6 7],squeeze(p(6:7)),'k.','markersize',25)
axis([0 14 -0.1 1])


% ATOMOXEINE

ax{3} = subplot(2,2,3); hold on
[~,p]=ttest(M_rest(:,2,1,:),M_rest(:,1,1,:))
plot(1:13,squeeze(p),'linewidth',3,'color','r')
line([0 13],[0.05 0.05],'linestyle',':','color',[0.7 0.7 0.7])
axis([0 14 -0.1 1])

[~,p]=ttest(M_rest(:,2,2,:),M_rest(:,1,2,:))
plot(1:13,squeeze(p),'linewidth',3,'color',[0.9 0.8 0.4])
line([0 13],[0.05 0.05],'linestyle',':','color',[0.7 0.7 0.7])
plot([6],squeeze(p(6)),'k.','markersize',25)
axis([0 14 -0.1 1])

ax{4} = subplot(2,2,4); hold on
[~,p]=ttest(M_rest(:,2,2,:)-M_rest(:,1,2,:),M_rest(:,2,1,:)-M_rest(:,1,1,:))
plot(1:13,squeeze(p),'linewidth',3,'color','k')
line([0 13],[0.05 0.05],'linestyle',':','color',[0.7 0.7 0.7])
plot([5 6 7],squeeze(p(5:7)),'k.','markersize',25)
axis([0 14 -0.1 1])

for iax = 1 : 4
    
  set(ax{iax},'XTick',1:2:13,'XTickLabel',[foi_range(1:2:13)])
  
  xlabel(ax{iax},'Frequency [Hz]')
  set(ax{iax},'YTick',[0 0.5 1],'YTickLabel',[0 0.5 1])
  ylabel(ax{iax},'P-Value (Uncorrected)')
  tp_editplots(ax{iax})
end

%% BARS FOR FC (HISTOGRAMS)

figure; set(gcf,'color','w'); tp_editplots; hold on
[n,k]=hist(fc(:,2,1,6)-fc(:,1,1,6),50);
bar(k,100*n./4005,'edgecolor','w','facecolor','b')
% plot(k,n,'r','linewidth',3)
[n,k]=hist(fc(:,2,2,6)-fc(:,1,2,6),50);
bar(k,100*n./4005,'edgecolor','w','facecolor','r')
axis([-0.01 0.01 -1 8]); axis square;
xlabel('Difference in FC'); ylabel('Number of connections')

figure; set(gcf,'color','w'); tp_editplots; hold on
[n,k]=hist(fc(:,3,1,7)-fc(:,1,1,7),50);
bar(k,100*n./4005,'edgecolor','w','facecolor','b')
% plot(k,n,'r','linewidth',3)
[n,k]=hist(fc(:,3,2,7)-fc(:,1,2,7),50);
bar(k,100*n./4005,'edgecolor','w','facecolor','r')
axis([-0.01 0.01 -1 8]); axis square;
xlabel('Difference in FC'); ylabel('Number of connections')

figure; set(gcf,'color','w'); tp_editplots; hold on
[n,k]=hist(fc(:,2,2,7),50);
bar(k,100*n./4005,'edgecolor','w','facecolor','b','facealpha',0.5)
% plot(k,n,'r','linewidth',3)
[n,k]=hist(fc(:,1,2,7),50);
bar(k,100*n./4005,'edgecolor','w','facecolor','r','facealpha',0.5)
axis([0 0.03 -1 8]); axis square;
xlabel('Difference in FC'); ylabel('Number of connections')

%%

d = fc(:,2,2,6)-fc(:,1,2,6);

[i,j]=sort(d,'descend')

k = j(1:100);
j = j(end:-1:end-100);

d1 = [mean(d(k)) mean(d(j))]
d2 = [mean(fc(k,2,1,6)-fc(k,1,1,6)) mean(fc(j,2,1,6)-fc(j,1,1,6))]

d = fc(:,3,1,6)-fc(:,1,1,6);

[i,j]=sort(d,'descend')

k = j(1:100);
j = j(end:-1:end-100);

d1 = [mean(d(k)) mean(d(j))]
d2 = [mean(fc(k,3,2,6)-fc(k,1,2,6)) mean(fc(j,3,2,6)-fc(j,1,2,6))]

% strongest vs. weakest
d = fc(:,1,1,6);
[i,j]=sort(d,'descend')
k = j(1:100);
j = j(end:-1:end-100);

d1 = mean(fc(k,2,2,6)-fc(k,1,2,6));
d2 = mean(fc(j,2,2,6)-fc(j,1,2,6));


