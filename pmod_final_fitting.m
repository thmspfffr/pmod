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

clear outp

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
fc_rest_indiv = reshape(squeeze(cleandat(:,:,:,1,1,6)),[90*90 28]);
fc_rest_indiv = fc_rest_indiv(mask,:);
fc_task_indiv = reshape(squeeze(cleandat(:,:,:,1,2,6)),[90*90 28]);
fc_task_indiv = fc_task_indiv(mask,:);
if ~exist(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_all_v%d.mat',v_sim))

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
        [outp.r_env_task_corr(iies,iiis,iG,igain), outp.p_env_task_corr(iies,iiis,iG,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_task(mask));

        [outp.r_rest_corr_avg(iies,iiis,iG,igain), outp.p_rest_corr_avg(iies,iiis,iG,igain)]=corr(nanmean(outp.fc_sim_tmp)',nanmean(fc_rest)');
        
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
  save(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_all_v%d.mat',v_sim),'outp')
% % % 
else
  load(sprintf('~/pmod/proc/pmod_wc_wholebrain_final_all_v%d.mat',v_sim))
end

error('!')
osc = osc1(:,:,4,1);

%% DETERMINE GLOBAL COUPLING PARAMETER FIRST

figure; set (gcf,'color','w');
subplot(2,2,1)
par1 = flipud(squeeze(mean(mean(outp.r_env_rest_corr(:,:,:,:))))');
imagesc(par1);
colormap(plasma);
tp_editplots
set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',num2cell(Gg))
set(gca,'YTick',1:5,'YTickLabel',num2cell(Gains(end:-1:1)));
xlabel('Global coupling'); ylabel('Gain')
cb = colorbar; set(cb,'Limits',[0 0.05],'Ticks',[0 0.05]); axis square
axis(gca,'square')

subplot(2,2,2)
par2 = flipud(squeeze(mean(mean(outp.fc_sim_env_mean(:,:,:,:))))');
imagesc(par2);
colormap(plasma);
tp_editplots
set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',num2cell(Gg))
set(gca,'YTick',1:5,'YTickLabel',num2cell(Gains(end:-1:1)));
xlabel('Global coupling'); %ylabel('Gain')
cb = colorbar; set(cb,'Limits',[0 0.18],'Ticks',[0 0.18]); axis square
axis(gca,'square')

subplot(2,2,3)
mask = par1 > 0.05 & par2 > 0.15;
imagesc(mask)
colormap(plasma);
tp_editplots
set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',num2cell(Gg))
set(gca,'YTick',1:5,'YTickLabel',num2cell(Gains(end:-1:1)));
xlabel('Global coupling'); %ylabel('Gain')
cb = colorbar; set(cb,'Limits',[0 0.18],'Ticks',[0 0.18]); axis square

axis(gca,'square')

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_final_fitG_v%d.pdf',v_sim))


%%
  load(sprintf('~/pmod/proc/pmod_wc_wholebrain_detosc_all_v%d.mat',v_sim))

osc = osc1(:,:,4);
oscthres = 0;
nancol = [0.97 0.97 0.97];

cmap = cbrewer('div', 'RdBu', 100,'pchip');
cmap = cmap(end:-1:1,:);

clear par ax
figure; set(gcf,'color','w')
% clear 

igain = 1;
iG = 4;
% plot peak freq model
ax{1} = subplot(2,2,1); hold on
par = squeeze(outp.peakfreq(:,:,iG,igain));
par(osc>oscthres)=nan;
imagescnan(par,[2 20],'NanColor',nancol)
title('Peak freq: Difference');

% plot peak freq model
ax{2} = subplot(2,2,2); hold on
par = squeeze(outp.peakfreq(:,:,iG,igain)-peakfreq_rest);
par(abs(par)<3)=1; par(abs(par)>3)=0; 
par(osc>oscthres)=0;
m1 = par>0;
par = double(par); par(par<1)=NaN;
imagescnan(par,[-1 1],'NanColor',nancol)
title('Peak freq: Masked');

% plot peak freq model
ax{3} = subplot(2,2,3); hold on
par = -log10(squeeze(outp.p_env_rest_corr(:,:,iG,igain)));
par(osc>oscthres)=nan;
imagescnan(par,[0 3],'NanColor',nancol)
title('r(FC_{sim},FC_{MEG})');

% plot peak freq model
fdr_p = fdr1(reshape(squeeze(outp.p_env_rest_corr(:,:,iG,igain)),[41*31 1]),0.05);
ax{4} = subplot(2,2,4); hold on
par = -log10(squeeze(outp.p_env_rest_corr(:,:,iG,igain)));
par(par<-log(fdr_p))=0; par(par>=-log(fdr_p))=1;
% par=par.*squeeze(outp.r_env_rest_corr(:,:,iG,igain))>0.15;
par(osc>oscthres)=0;
m2 = par>0;
par = double(par); par(par<1)=NaN;
imagescnan(par,[-1 1],'NanColor',nancol)
title('Peak freq: Masked');

[i,j]= find(m2);
idx = round([mean(i) mean(j)]);

for iax = 1 : 4
 scatter(ax{4},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
%   scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','r','markeredgecolor','k')
  if iax == 1
    colormap(ax{iax},plasma)
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  elseif iax == 2
    colormap(ax{iax},redblue)
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  elseif iax == 3
    colormap(ax{iax},plasma)
    xlabel(ax{iax},'Inhibitory input')
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  elseif iax ==4
    colormap(ax{iax},redblue)
    xlabel(ax{iax},'Inhibitory input')
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  end
  tp_editplots(ax{iax})
  
  c = colorbar(ax{iax}); axis(ax{iax},'tight')
  c.Ticks = [min(c.Limits) max(c.Limits)];
  if iax == 3
    c.TickLabels = {'0'; sprintf('%.3f',10^-max(c.Limits))};
  end
   axis(ax{iax},[1 length(Iis) 1 length(Ies) ])
  c.Ticks = [min(c.Ticks) max(c.Ticks)];
   axis(ax{iax},'square')
  
end

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_final_fitting_fc_v%d.pdf',v_sim))

%% FIT TASK PARAMETERS (MEAN - see below for indiv subjects)

mask = logical(tril(ones(90,90),-1));
fc_tvr_sim = 100*(outp.fc_sim_mean(:,:,4,1)-outp.fc_sim_mean(idx(1),idx(2),4,1))./outp.fc_sim_mean(idx(1),idx(2),4,1)
fc_tvr_exp = 100*(mean(fc_task(mask))-mean(fc_rest(mask)))./mean(fc_rest(mask));

% plot peak freq model
ax{4} = subplot(2,2,4); hold on
par = fc_tvr_sim;
par(par<(fc_tvr_exp-10))=0; par(par>0)=0; par(par<(fc_tvr_exp+10))=1; 
par(1:idx(1),1:idx(2))=nan
% par=par.*squeeze(outp.r_env_rest_corr(:,:,iG,igain))>0.15;
par(osc>oscthres)=0;
m2 = par>0;
par = double(par); par(par<1)=NaN;
imagescnan(par,[-1 1],'NanColor',nancol)
title('FC decr. Task (+/- 10%)');
colormap(redblue)

scatter(gca,idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')

[idx_task(:,1),idx_task(:,2)]=find(par>0)

scatter(gca,idx_task(:,2),idx_task(:,1),20,'markerfacecolor','r','markeredgecolor','k')

iax = 4
colormap(ax{iax},redblue)
ylabel(ax{iax},'Excitatory input');
set(ax{iax},'YTick',1:50:length(Ies ),'YTickLabels',num2cell(Ies(1:50:end)))
set(ax{iax},'XTick',1:100:length(Iis),'XTickLabels',num2cell(Iis(1:100:end)))
xlabel(ax{iax},'Inhibitory input'); tp_editplots;
axis(ax{iax},[1 length(Iis) 1 length(Ies) ])
%   c.Ticks = [min(c.Ticks) max(c.Ticks)];
axis(ax{iax},'square')

clear par
par.rest = [Ies(idx(1)) Iis(idx(2))];
par.task = [Ies(idx_task(1,1)) Iis(idx_task(1,2))];
par.task_alt = [Ies(idx_task(2:end,1)); Iis(idx_task(2:end,2))]';
par.descr = 'First entry: Excitation (Ies), Second entry: Inhibition (Iis)';

save(sprintf('~/pmod/proc/pmod_final_fitting_fits_v%d.mat',v_sim),'par')


%%
[i,k]=find(m1&m2)
idx = [round(mean(i)) round(mean(k))];

figure; set(gcf,'color','w')
ax{1} = subplot(2,2,1); hold on
par = double(m1&m2);
par(osc>oscthres)=nan;
imagescnan(par,[-1 1])
title('Peak freq: Masked');
colormap(gca,redblue)

d_frest_ftask = peakfreq_task-peakfreq_rest;

% plot peak freq model
ax{2} = subplot(2,2,2); hold on
par = (squeeze(outp.peakfreq(:,:,4,1))-outp.peakfreq(idx(1),idx(2),4,1))-d_frest_ftask;
par(abs(par)<2)=1; par(abs(par)>=2)=0;
par(osc>oscthres)=nan; 
m1 = par>0;
imagescnan(par,[-1 1])
title('Peak freq: Task vs rest');

% plot peak freq model
ax{3} = subplot(2,2,3); hold on
par = double(m1&m2);
par(osc>oscthres)=nan;

imagescnan(par,[-1 1])
title('Peak freq: Difference');
colormap(redblue)


[i,k]=find(m1&m2)
idx2 = [round(mean(i)) round(mean(k))];

ax{4} = subplot(2,2,4); hold on
par =zeros(size((par)));
% par(osc1>0.5)=nan;
imagescnan(par,[-1 1]);
title('Peak freq: Difference');


for iax = 1 : 4
  scatter(ax{iax},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
  if iax == 1
    colormap(ax{iax},redblue)
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  elseif iax == 2
    colormap(ax{iax},redblue)
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  elseif iax == 3
    colormap(ax{iax},redblue)
    xlabel(ax{iax},'Inhibitory input')
    ylabel(ax{iax},'Excitatory input');
    set(ax{iax},'YTick',1:5:length(Ies ),'YTickLabels',num2cell(Ies(1:5:end)))
    set(ax{iax},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
  elseif iax == 4
    colormap(ax{iax},redblue)
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

save(sprintf('~/pmod/proc/pmod_final_fitting_fits_v%d.mat',v_sim),'par')



%% INDIV SUBJ

% osc=zeros(301,401)

h=figure; set(gcf,'color','w')

width = 4;
bif_mask = zeros(31,41);
for i = 1 : 31
  if isempty(find(osc(i,:)>0,1,'last'))
    idx=0;
  else
    idx = find(osc(i,:)>0,1,'last')
    if idx-width < 1

      bif_mask(i,idx-(width+(idx-(width+1))):idx+width)=1;
    else
      bif_mask(i,idx-width:idx+width)=1;
    end
  end
end
clear idx
ax{1}=subplot(6,5,1); hold on

imagescnan(bif_mask,[-1 1])

colormap(redblue)
set(gca,'YTick',1:10:length(Ies ),'YTickLabels',num2cell(Ies(1:10:end)))
set(gca,'YTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
tp_editplots

for isubj = 1 : 28
  % plot lambda
  ax{1} = subplot(6,5,isubj+1); hold on
  par = squeeze(outp.p_env_rest_indiv_corr(isubj,:,:,1,1))<0.05; par=double(par);
  
  par(osc>0)=0;
  par(~bif_mask)=0
  
  bw = bwlabel(par,8);
  bw(osc>0)=0;
  for i = 1 : max(bw(:))
    if sum(bw(:)==i) <= 3
      bw(bw==i)=0;
    end
  end

  par = bw;
  par(osc>0)=nan;

  imagescnan(par,[-1 1])
  colormap(redblue)

  set(gca,'YTick',1:100:length(Ies ),'YTickLabels',num2cell(Ies(1:100:end)))
  set(gca,'XTick',1:100:length(Iis),'XTickLabels',num2cell(Iis(1:100:end)))

  [i,k]=find(par>0);
  idx(isubj,:) = [round(mean(i)) round(mean(k))];
  scatter(gca,idx(isubj,2),idx(isubj,1),20,'markerfacecolor','w','markeredgecolor','k');
  tp_editplots
end

set(h,'Position',[50 50 1200 900])
set(h,'Renderer','Painters')
print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_fits_indiv.pdf'))

save(sprintf('~/pmod/proc/pmod_final_fitting_indivfits_v%d.mat',v_sim),'idx')

close all
%%
figure; set(gcf,'color','white')

ax{1} = subplot(2,2,1); hold on

imagescnan(bif_mask,[-1 1])

colormap(redblue)
set(gca,'YTick',1:10:length(Ies ),'YTickLabels',num2cell(Ies(1:10:end)))
set(gca,'YTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
tp_editplots
title('Bifurcation Mask')
set(ax{1},'YTick',1:10:length(Ies ),'YTickLabels',num2cell(Ies(1:10:end)))
set(ax{1},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
tp_editplots
xlabel(ax{1},'Inhibitory input')
ylabel(ax{1},'Excitatory input');
axis(ax{1},'square')
colorbar
ax{1} = subplot(2,2,4); hold on
scatter(gca,idx(:,2),idx(:,1),10,'markerfacecolor','w','markeredgecolor','k');
axis(gca,[1 41 1 31])

set(ax{1},'YTick',1:10:length(Ies ),'YTickLabels',num2cell(Ies(1:10:end)))
set(ax{1},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
tp_editplots
xlabel(ax{1},'Inhibitory input')
ylabel(ax{1},'Excitatory input');
axis(ax{1},'square')
title('Indiv subjects')
colorbar
print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_fits_indiv_sum.pdf'))


%% GET INDIV TASK PARAMETERS
indiv_change_prc = 100*(mean(fc_task_indiv,1)-mean(fc_rest_indiv,1))./mean(fc_rest_indiv,1);

for isubj = 1:28
  if isnan(idx(isubj,1))
    continue
  end
  fc_tvr_sim = 100*(outp.fc_sim_mean(:,:,4,1)-outp.fc_sim_mean(idx(isubj,1),idx(isubj,2),4,1))./outp.fc_sim_mean(idx(isubj,1),idx(isubj,2),4,1)

  d = indiv_change_prc(isubj)-fc_tvr_sim;
  d(osc>0) = nan;



thresh = 1;
d(~bif_mask) = Inf;
d(1:idx(isubj,1),1:idx(isubj,2)) = Inf;

while 1
  clust = reshape(bwlabel(abs(d)<thresh,8),[31*41 1]);
  if max(clust)<1
    thresh = thresh+1;
    continue
  else
    clear cnt
    for iclust = 1 : max(clust)
      cnt(iclust) = sum(clust==iclust);
    end
    if any(cnt>1)
      break
    else
      thresh = thresh+1;
      continue
    end
  end
end

[i,j] = find(bwlabel(abs(d)<thresh)==find(cnt>1,1,'first'))
idx_task(isubj,:) = round([mean(i) mean(j)]);

end


figure; set(gcf,'color','w')

ax{1}= subplot(2,2,1); hold on

imagesc(zeros(31,41))
scatter(ax{1},idx_task(:,2),idx_task(:,1),20,'markerfacecolor','y','markeredgecolor','k');
set(ax{1},'YTick',1:10:length(Ies),'YTickLabels',num2cell(Ies(1:10:end)))
set(ax{1},'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
xlabel(ax{1},'Inhibitory input'); ylabel(ax{1},'Excitatory input');
axis(ax{1},'square')
axis([1 41 1 31])
tp_editplots; colormap(redblue); colorbar
print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_fits_indiv_task_sum.pdf'))

load(sprintf('~/pmod/proc/pmod_final_fitting_indivfits_v%d.mat',v_sim))

clear par

indiv_idx.rest = idx;
indiv_idx.task = idx_task;

save(sprintf('~/pmod/proc/pmod_final_fitting_indivfits_taskandrest_v%d.mat',v_sim),'indiv_idx')

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

