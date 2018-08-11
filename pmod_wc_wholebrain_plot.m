%% FITTING
% pmod_wc_wholebrain_plot.m
Ies         = -4:0.1:-1;
Iis         = -5:0.1:-1;
Gg          = 0.62;
Gains       = -0.5:0.25:0.5;
% connectivity, AAL
v_conn =  1;
% simulations
v_sim = 2;
% dfa, aal, lcmv
v_dfa = 2;
% -------------------

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

clear dfa r dfa_r dist_fc fc_sim autocorr dfa_sim dfa_env_sim r*

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

if ~exist(sprintf('~/pmod/proc/pmod_wholebrain_rest_all_v%d.mat',v_sim))

for iies = 1 : length(Ies)
  iies
  for iiis = 1 : length(Iis)
    for iG = Gs
      for igain = 1:5
        %         igain
        load(sprintf('~/pmod/proc/pmod_wc_wholebrain_rest_Ie%d_Ii%d_G%d_gain%d_v%d.mat',iies,iiis,iG,igain,v_sim))

        if round(Ies(iies)*100)/100 == -2.8 && round( Iis(iiis)*100)/100 == -3.4
          disp('save stuff')
          FFCC = out.FC;
        end
        
        % Time scales
        outp.lambda(:,iies,iiis,iG,igain)      = mean(out.lambda,2);
        outp.lambda_env(:,iies,iiis,iG,igain)  = mean(out.lambda_env,2);
        % DFA
        outp.dfa_sim(:,iies,iiis,iG,igain)     = squeeze(mean(out.dfa,1));
        outp.dfa_env_sim(:,iies,iiis,iG,igain) = squeeze(mean(out.dfa_env,1));
        
        [~,peak_idx]=max(smooth(mean(out.PSD(out.f>3,:),2),20));
        outp.peakfreq(iies,iiis,iG,igain) = out.f(peak_idx+find(out.f<3,1,'last'));
        
        % get peak frequency from MEG data
        outp.peakfreq_diff_res(iies,iiis,iG,igain) = outp.peakfreq(iies,iiis,iG,igain)-peakfreq_rest;
        outp.peakfreq_diff_tsk(iies,iiis,iG,igain) = outp.peakfreq(iies,iiis,iG,igain)-peakfreq_task;
        
        pars.dim = 2;
        
        outp.fc_sim_tmp = tp_match_aal(pars,out.FC,3);     
        outp.fc_sim_env = tp_match_aal(pars,out.FC_env,3);
        
        outp.fc_sim_var(:,iies,iiis,iG,igain) = std(nanmean(outp.fc_sim_tmp)./max(nanmean(outp.fc_sim_tmp)));
        
        [outp.r_rest_corr(iies,iiis,iG,igain), outp.p_rest_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_tmp(mask),fc_rest(mask));
        [outp.r_rest_corr_avg(iies,iiis,iG,igain), outp.p_rest_corr_avg(iies,iiis,iG,igain)]=corr(nanmean(outp.fc_sim_tmp)',nanmean(fc_rest)');
        
        outp.r_rest_corr_unc(iies,iiis,iG,igain) = dot(outp.fc_sim_tmp(mask),fc_rest(mask)) / sqrt(dot(outp.fc_sim_tmp(mask),outp.fc_sim_tmp(mask)) * dot(fc_rest(mask),fc_rest(mask)));
        %
        pars.dim = 1;
        
        [outp.dfa_r_rest(iies,iiis,iG,igain), outp.dfa_p_rest(iies,iiis,iG,igain)] = corr(dfa_emp_rest(:),tp_match_aal(pars,repmat(outp.dfa_sim(:,iies,iiis,iG,igain),[1 90]),pars));
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
        
        outp.kuramoto_mean_diff (iies,iiis,iG,igain) = mean(out.KOPmean) - mean(kura_mean(:,1,1,1));
        outp.kuramoto_std_diff (iies,iiis,iG,igain)  = mean(out.KOPsd)- mean(kura_std(:,1,1,1));
        
        outp.psslp(:,iies,iiis,iG,igain)      = out.psslope;
        outp.psslp_env(:,iies,iiis,iG,igain)  = out.psslope_env;
        
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
% % 
else
  load(sprintf('~/pmod/proc/pmod_wholebrain_rest_all_v%d.mat',v_sim))
end

error('!')
%% PLOT BASIC PARAMETERS
% Lambda
% Peak frequency
% Correlation model with MEG: Lambda and FC
% ------------------------------------------
clear par

idx = [find(round(Ies*100)/100==-2.8) find(round(Iis*100)/100==-3.5000)];
idx2 = [find(round(Ies*100)/100==-1.8) find(round(Iis*100)/100==-2.4000)];

igain = 3;
G = 1;

figure; set(gcf,'color','w')

% plot lambda
ax{1} = subplot(2,2,1); hold on
par = squeeze(mean(1./outp.lambda(:,:,:,G,igain)));
par(osc1>0.5)=nan;
imagescnan(par,[4 18])
title('\lambda (model)');

% plot correlation FC model / MEG
ax{2} = subplot(2,2,2); hold on
par = squeeze(outp.r_rest_corr(:,:,G,igain));
par(osc1>0.5)=nan;
imagescnan(par,[0 0.3])
title('Corr(FC_{sim}, FC_{exp})');

% plot correlation lambda model / MEG
ax{3} = subplot(2,2,3); hold on
par = squeeze(outp.lambda_r_rest(:,:,G,igain));
par(osc1>0.5)=nan;
imagescnan(par,[0 0.10])
title('Corr(\lambda_{sim},\lambda_{exp})');

% plot peak freq model
ax{4} = subplot(2,2,4); hold on
par = squeeze(outp.peakfreq(:,:,G,igain));
par(osc1>0.5)=nan;
imagescnan(par,[2 20])
title('Peak frequency');

for iax = 1 : length(ax)
  
  scatter(ax{iax},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
  scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
  ylabel(ax{iax},'Excitatory input'); xlabel(ax{iax},'Inhibitory input')
  set(ax{iax},'XTick',1:5:length(Iis),'XTickLabels',Iis(1:5:end))
  set(ax{iax},'YTick',1:5:length(Ies),'YTickLabels',Ies(1:5:end))
  tp_editplots(ax{iax})
  colormap(plasma)
  
end

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_lambda_gain%d_G%d_v%d.pdf',igain,G,v_sim))


%% PLOT GAIN VS NO GAIN
% increase in gain vs. baseline
clear par

idx = [find(round(Ies*100)/100==-2.8) find(round(Iis*100)/100==-3.5000)];
idx2 = [find(round(Ies*100)/100==-1.8) find(round(Iis*100)/100==-2.4000)];

igain = 4;
G = 1;

cmap = cbrewer('div', 'RdBu', 100,'pchip');
cmap = cmap(end:-1:1,:);

figure; set(gcf,'color','w')

% plot lambda
ax{1} = subplot(2,2,1); hold on
par = squeeze(outp.fc_sim_mean(:,:,G,3));
par(osc1>0.5)=nan;
imagescnan(par,[0 0.02])
title('Mean FC_{sim}');

% plot correlation FC model / MEG
ax{2} = subplot(2,2,2); hold on
par = squeeze(abs(outp.fc_sim_mean(:,:,G,igain)))-abs(squeeze(outp.fc_sim_mean(:,:,G,3)));
par(osc1>0.5)=nan;
imagescnan(par,[-0.02 0.02])
title('FC_{sim}: Gain vs Baseline');

% plot correlation lambda model / MEG
ax{3} = subplot(2,2,3); hold on
par = squeeze(mean(1./outp.lambda(:,:,:,G,igain)))-squeeze(mean(1./outp.lambda(:,:,:,G,3)));
par(osc1>0.5)=nan;
imagescnan(par,[-10 10])
title('\lambda: Gain vs Baseline');

% plot peak freq model
ax{4} = subplot(2,2,4); hold on
par = squeeze(outp.peakfreq(:,:,G,igain))-squeeze(outp.peakfreq(:,:,G,3));
par(osc1>0.5)=nan;
imagescnan(par,[-10 10])
title('Peak frequency');

for iax = 1 : length(ax)
  
  scatter(ax{iax},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
  scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
  ylabel(ax{iax},'Excitatory input'); xlabel(ax{iax},'Inhibitory input')
  set(ax{iax},'XTick',1:5:length(Iis),'XTickLabels',Iis(1:5:end))
  set(ax{iax},'YTick',1:5:length(Ies),'YTickLabels',Ies(1:5:end))
  tp_editplots(ax{iax})
  colormap(ax{1},plasma)
  colormap(ax{iax},cmap)
  
end

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_gainvsbaseline_gain%d_G%d_v%d.pdf',igain,G,v_sim))


%%
clear par
figure; set(gcf,'color','w')
clear ax

igain = 3;
% plot peak freq model
ax{1} = subplot(2,2,1); hold on
par = squeeze(outp.degree(:,:,G,igain));
par(osc1>0.5)=nan;
imagescnan(par,[0 1])
title('Degree');

% plot peak freq model
ax{2} = subplot(2,2,2); hold on
par = squeeze(outp.degree(:,:,G,4))-squeeze(outp.degree(:,:,G,3));
par(osc1>0.5)=nan;
imagescnan(par,[-0.25 0.25])
title('Difference Degree');

for iax = 1 : length(ax)
  
  scatter(ax{iax},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
  scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
  ylabel(ax{iax},'Excitatory input'); xlabel(ax{iax},'Inhibitory input')
  set(ax{iax},'XTick',1:5:length(Iis),'XTickLabels',Iis(1:5:end))
  set(ax{iax},'YTick',1:5:length(Ies),'YTickLabels',Ies(1:5:end))
  tp_editplots(ax{iax})
  colormap(ax{1},plasma)
  colormap(ax{iax},cmap)
  axis(ax{iax}, 'tight')
end

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_degree_gain%d_G%d_v%d.pdf',igain,G,v_sim))

%% KURAMOTO

clear par ax

figure; set(gcf,'color','w')

igain = 3; G = 1;
% Kuramoto mean
ax{1} = subplot(2,2,1); hold on
par = squeeze(outp.kuramoto_mean(:,:,G,igain));
par(osc1>0.5)=nan;
imagescnan(par,[0 0.25])
title('Kuramoto parameter');

% Kuramoto std - metastability
ax{2} = subplot(2,2,2); hold on
par = squeeze(outp.kuramoto_std(:,:,G,igain));
par(osc1>0.5)=nan;
imagescnan(par,[0 0.0002])
title('Metastability');

% Kuramoto mean: contrast
ax{3} = subplot(2,2,3); hold on
par = squeeze(outp.kuramoto_mean(:,:,G,4)-outp.kuramoto_mean(:,:,G,igain));
par(osc1>0.5)=nan;
imagescnan(par,[-0.075 0.075])
title('Difference Kuramoto');

% Kuramoto std - metastability
ax{4} = subplot(2,2,4); hold on
par = squeeze(outp.kuramoto_std(:,:,G,4)-outp.kuramoto_std(:,:,G,igain));
par(osc1>0.5)=nan;
imagescnan(par,[-0.0001 0.0001])
title('Difference Metastability');


for iax = 1 : length(ax)
  
  scatter(ax{iax},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
  scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
  ylabel(ax{iax},'Excitatory input'); xlabel(ax{iax},'Inhibitory input')
  set(ax{iax},'XTick',1:5:length(Iis),'XTickLabels',Iis(1:5:end))
  set(ax{iax},'YTick',1:5:length(Ies),'YTickLabels',Ies(1:5:end))
  tp_editplots(ax{iax})
  colormap(ax{1},plasma)
  colormap(ax{2},plasma)
  colormap(ax{iax},cmap)
  
%   axis(ax{iax}, 'tight')
end

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_kuramoto_gain%d_G%d_v%d.pdf',igain,G,v_sim))


%%
drug_gain = 4;

figure; set(gcf,'color','w')

subplot(1,2,1);

fc = [squeeze(outp.fc_sim_mean(idx(1),idx(2),G,3)),squeeze(outp.fc_sim_mean(idx(1),idx(2),G,drug_gain))];
fc = [fc squeeze(outp.fc_sim_mean(idx2(1),idx2(2),G,3)),squeeze(outp.fc_sim_mean(idx2(1),idx2(2),G,drug_gain))];

bar([1 2 3 4],fc);
axis([0 5 0 0.008])
tp_editplots; ylabel('Mean correlation'); axis square

subplot(1,2,2);
deg = [squeeze(outp.degree(idx(1),idx(2),G,3)),squeeze(outp.degree(idx(1),idx(2),G,drug_gain))];
deg = [deg squeeze(outp.degree(idx2(1),idx2(2),G,3)),squeeze(outp.degree(idx2(1),idx2(2),G,drug_gain))];

bar([1 2 3 4 ],[100.*deg]);
axis([0 5 25 70])
tp_editplots; ylabel('Degree [%]'); axis square

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_degree_gain%d_G%d_v%d.pdf',drug_gain,G,v_sim))

%% TESTS

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

