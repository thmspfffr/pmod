%% FITTING
% pmod_wc_ei_plot.m
%-------------------------------------------------------------------------
% VERSION 2: After meeting with tobi, 24-08-2018
% This is based on the fine grained simulations and indivudual fits
% obtained in pmod_wc_highres_fitting.m
%-------------------------------------------------------------------------
load(sprintf('~/pmod/proc/pmod_wc_highres_indivfits_taskandrest_v%d.mat',3))
% parameters are also saved in home folder
v           = 2;
Ies         = indiv_idx.ie;
Iis         = indiv_idx.ii;
EIs         = 0.8:0.02:1.2;
Gains       = -0.2:0.02:0.2;
nTrials     = 1;
tmax        = 6500; % in units of tauE
Gg          = 0.60;
wins        = [2 20];
%-------------------------------------------------------------------------
% VERSION 3: more subtle changes
%-------------------------------------------------------------------------
load(sprintf('~/pmod/proc/pmod_wc_highres_indivfits_taskandrest_v%d.mat',3))
v           = 3;
Ies         = indiv_idx.ie;
Iis         = indiv_idx.ii;
EIs         = 0.95:0.005:1.05;
Gains       = -0.1:0.01:0.1;
nTrials     = 1;
tmax        = 6500; % in units of tauE
Gg          = 0.60;
wins        = [2 20];
%-------------------------------------------------------------------------
% VERSION 4: CHANGES IN E/I separately
%-------------------------------------------------------------------------
load(sprintf('~/pmod/proc/pmod_wc_highres_indivfits_taskandrest_v%d.mat',3))
v           = 4;
Ies         = [round(mean(indiv_idx.ie(1:28)*10))/10 round(mean(indiv_idx.ie(29:56)*10))/10];
Iis         = [round(mean(indiv_idx.ii(1:28)*10))/10 round(mean(indiv_idx.ii(29:56)*10))/10];
wEE_mod     = 0.8:0.01:1.2;
wIE_mod     = 0.8:0.01:1.2;
nTrials     = 1;
tmax        = 6500; % in units of tauE
Gg          = 0.60;
wins        = [2 20];
Gains       = 0;
%-------------------------------------------------------------------------

v_sim = v;
% connectivity, AAL
v_conn =  1;
% simulations
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
pars.grid = 'medium';-0.1:0.01:0.1
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

pars          = [];
pars.dim      = 2;
pars.grid     = 'medium';
pars.N        = 90;
pars.transfer = 'to_bcn';
fc_rest     =  tp_match_aal(pars,squeeze(nanmean(cleandat(:,:,:,1,1,6),3)));
fc_task     =  tp_match_aal(pars,squeeze(nanmean(cleandat(:,:,:,1,2,6),3)));


% if ~exist(sprintf('~/pmod/proc/pmod_wc_ei_poi_all_v%d.mat',v_sim))
for iies = 1: length(Ies)
  iies
  for EE = 1 : length(wEE_mod)
    for IE = 1 : length(wIE_mod)
      for igain = 1 : length(Gains)
        
        fn = sprintf('pmod_wc_ei_poi_Ie%d_EE%d_IE%d_gain%d_v%d',iies,EE,IE,igain,v);
       
        
        load(['~/pmod/proc/ei/' fn '.mat'])

        outp.peakfreq(iies,EE,IE,igain) = out.peakfreq;
                
        outp.fc_sim_tmp = out.FC;
        outp.fc_sim_env_tmp = out.FC_env;
        
        [outp.r_rest_corr(iies,EE,IE,igain), outp.p_rest_corr(iies,EE,IE,igain)]=corr(outp.fc_sim_tmp(mask),fc_rest(mask));
        [outp.r_env_rest_corr(iies,EE,IE,igain), outp.p_env_rest_corr(iies,EE,IE,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_rest(mask));
        [outp.r_env_task_corr(iies,EE,IE,igain), outp.p_env_task_corr(iies,EE,IE,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_task(mask));
                
        pars.dim = 1;

        outp.dist_fc_rest (iies,EE,IE,igain)  = mean(outp.fc_sim_tmp(mask))-mean(fc_rest(mask));
        outp.dist_fc_task (iies,EE,IE,igain)  = mean(outp.fc_sim_tmp(mask))-mean(fc_task(mask));
        
        outp.fc_sim_mean(iies,EE,IE,igain) = mean(outp.fc_sim_tmp(mask));
        outp.fc_sim_env_mean(iies,EE,IE,igain) = mean(outp.fc_sim_env_tmp(mask));

        outp.EIratio(iies,EE,IE,igain) = mean(out.EIratio);
        
        fclose all;
        %
      end
    end
  end
end
  save(sprintf('~/pmod/proc/ei/pmod_wc_ei_poi_all_v%d.mat',v_sim),'outp')
  % % %
% else
%   load(sprintf('~/pmod/proc/pmod_wc_ei_poi_all_v%d.mat',v_sim))
% end

error('!')
%%
figure; set(gcf,'color','w');  hold on
ac{1}=subplot(1,2,1);hold on

imagesc(squeeze(outp.EIratio(1,:,:,11)),[0.75 1.25])
xlabel(gca,'wIE'); ylabel(gca,'wEE');
set(ac{1},'YTick',1:10:41,'YTickLabels',num2cell(wEE_mod(1:10:41)))
set(ac{1},'XTick',1:10:41,'XTickLabels',num2cell(wIE_mod(1:10:41)))
tp_editplots
scatter(21,21,50,'markerfacecolor','w','markeredgecolor','w')
axis([1 41 1 41]); axis square
colorbar

ac{2}=subplot(1,2,2);hold on
imagesc(squeeze(outp.EIratio(1,:,:,11))-squeeze(outp.EIratio(1,21,21,11)),[-0.2 0.2])
xlabel(gca,'wIE'); ylabel(gca,'wEE');
set(ac{2},'YTick',1:10:41,'YTickLabels',num2cell(wEE_mod(1:10:41)))
set(ac{2},'XTick',1:10:41,'XTickLabels',num2cell(wIE_mod(1:10:41)))
tp_editplots
scatter(21,21,50,'markerfacecolor','w','markeredgecolor','w')
axis([1 41 1 41]); axis square
colorbar

colormap(ac{1},plasma)
colormap(ac{2},redblue)

%%
figure; set(gcf,'color','w');  hold on
ac{1}=subplot(1,2,1);hold on

imagesc(squeeze(outp.fc_sim_env_mean(1,:,:,11)),[0 1])
xlabel(gca,'wIE'); ylabel(gca,'wEE');
set(ac{1},'YTick',1:10:41,'YTickLabels',num2cell(wEE_mod(1:10:41)))
set(ac{1},'XTick',1:10:41,'XTickLabels',num2cell(wIE_mod(1:10:41)))
tp_editplots
scatter(21,21,50,'markerfacecolor','w','markeredgecolor','w')
axis([1 41 1 41]); axis square
colorbar

ac{2}=subplot(1,2,2);hold on
imagesc(squeeze(outp.fc_sim_env_mean(1,:,:,11))-squeeze(outp.fc_sim_env_mean(1,21,21,11)),[-0.5 0.5])
xlabel(gca,'wIE'); ylabel(gca,'wEE');
set(ac{2},'YTick',1:10:41,'YTickLabels',num2cell(wEE_mod(1:10:41)))
set(ac{2},'XTick',1:10:41,'XTickLabels',num2cell(wIE_mod(1:10:41)))
tp_editplots
scatter(21,21,50,'markerfacecolor','w','markeredgecolor','w')
axis([1 41 1 41]); axis square
colorbar

colormap(ac{1},plasma)
colormap(ac{2},redblue)


%% PLOT GAIN VS NO GAIN
% increase in gain vs. baseline
clear par
% load(sprintf('~/pmod/proc/pmod_final_fitting_fits_v%d.mat',4))
load(sprintf('~/pmod/proc/pmod_final_fitting_indivfits_taskandrest_v%d.mat',4))

indiv_idx.task =indiv_idx.task(~isnan(indiv_idx.rest(:,1)),:);
indiv_idx.rest =indiv_idx.rest(~isnan(indiv_idx.rest(:,1)),:);

oscthres = 0
idx = [indiv_idx.rest(:,1) indiv_idx.rest(:,2)];
idx2 = [indiv_idx.task(:,1) indiv_idx.task(:,2)];
% [i,idx3(:,1)]=find(Ies==par.task_alt(:,1));
% [i,idx3(:,2)]=find(Iis==par.task_alt(:,2));

nancol = [0.95 0.95 0.95]
iei = 6;

% cmap = cbrewer('div', 'RdBu', 256,'pchip');
% cmap = cmap(end:-1:1,:);

figure; set(gcf,'color','w')

% plot lambda
ax{1} = subplot(2,2,1); hold on
par = squeeze(abs(outp.fc_sim_mean(:,:,iei)))-abs(squeeze(outp.fc_sim_mean(:,:,7)));
par(osc>oscthres)=nan;
imagescnan(par,[-0.02 0.02],'NanColor',nancol)
title('Contrast: FC_{FR}');

% plot correlation lambda model / MEG
ax{2} = subplot(2,2,2); hold on
par = squeeze(abs(outp.fc_sim_env_mean(:,:,iei)))-abs(squeeze(outp.fc_sim_env_mean(:,:,7)));
par(osc>oscthres)=nan;
imagescnan(par,[-0.02 0.02],'NanColor',nancol)
title('Contrast: FC_{env}');

ax{3} = subplot(2,2,3); hold on
par = squeeze(mean(1./outp.lambda(:,:,:,iei)))-squeeze(1./mean(outp.lambda(:,:,:,7)));
par(osc>oscthres)=nan;
imagescnan(par,[-5 5],'NanColor',nancol)
title('Contrast: Lambda_{FR}');

% plot peak freq model
ax{4} = subplot(2,2,4); hold on
par = squeeze(1./mean(outp.lambda_env(:,:,:,iei)))-squeeze(1./mean(outp.lambda_env(:,:,:,7)));
par(osc>oscthres)=nan;
imagescnan(par,[-30 30],'NanColor',nancol)
title('Contrast: Lambda_{Env}');

for iax = 1 : length(ax)
  scatter(ax{iax},idx(:,2),idx(:,1),20,'markerfacecolor','w','markeredgecolor','k')
  scatter(ax{iax},idx2(:,2),idx2(:,1),20,'markerfacecolor','y','markeredgecolor','k')
  c = colorbar(ax{iax});
  c.Ticks = [c.Limits];
  %   scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
  if iax == 1
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
    %     c.Position = [0.44 0.588 0.01 0.3310];
  elseif iax == 2
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
    %     c.Position = [0.88 0.588 0.01 0.3310];
  elseif iax == 3
    xlabel(ax{iax},'Inhibitory input')
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
    %     c.Position = [0.44 0.114 0.01 0.3310];
  elseif iax ==4
    xlabel(ax{iax},'Inhibitory input')
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
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


print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_ei_eivsbaseline_ei%d_v%d.pdf',iei,v_sim))
%%

% BAR PLOTS
tmp1 = squeeze(abs(outp.fc_sim_env_mean(1:28,:,:)))-squeeze(abs(outp.fc_sim_env_mean(1:28,11,11)));
tmp2 = squeeze(abs(outp.fc_sim_env_mean(29:56,:,:)))-squeeze(abs(outp.fc_sim_env_mean(29:56,11,11)));

diag_mask = eye(size(tmp1,1),size(tmp1,2));

tmp1 = squeeze(mean(tmp1));
tmp2 = squeeze(mean(tmp2));

par1 = tmp1(logical(diag_mask));
par2 = tmp2(logical(diag_mask));


figure;set(gcf,'color','w')
subplot(2,2,1); hold on
m = [mean(par1(:,iei)) mean(par2(:,iei))];
s = [std(par1(:,iei))/sqrt(length(par1)) std(par2(:,iei))/sqrt(length(par1))];
bar([1:2],[m])
line([1 1],[m(1)-s(1) m(1)+s(1)])
line([2 2],[m(2)-s(2) m(2)+s(2)])

tp_editplots; axis square
axis([0 3 -0.03 0.01])

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_ei_eivsbaseline_bar_ei%d_v%d.pdf',iei,v_sim))


figure;set(gcf,'color','w')
subplot(2,2,1); hold on
scatter(ones(size(par1,1),1)-(0.5-rand(22,1))/2,par1(:,iei),30,'markeredgecolor','k','markerfacecolor','w')
scatter(2*ones(size(par2,1),1)-(0.5-rand(22,1))/2,par2(:,iei),30,'markeredgecolor','k','markerfacecolor','y')


tp_editplots; axis square
axis([0 3 -0.02 0.02])

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_ei_eivsbaseline_scatter_ei%d_v%d.pdf',iei,v_sim))

%%



