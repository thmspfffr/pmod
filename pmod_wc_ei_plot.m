%% FITTING
% pmod_wc_ei_plot.m
%-------------------------------------------------------------------------
% VERSION 1: After meeting with tobi, 15-08-2018
%-------------------------------------------------------------------------
v           = 1;
Ies         = -4:0.1:-1;
Iis         = -5:0.1:-1;
EIs         = 0.7:0.05:1.3;
out.Gain = 0;
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

fc_rest     =  squeeze(nanmean(cleandat(:,:,:,1,1,6),3));
fc_task     =  squeeze(nanmean(cleandat(:,:,:,1,2,6),3));
fc_rest_var =  std(nanmean(squeeze(nanmean(cleandat(:,:,:,1,1,6),3)))./max(nanmean(squeeze(nanmean(cleandat(:,:,:,1,1,6),3)))));

if ~exist(sprintf('~/pmod/proc/pmod_wc_ei_all_v%d.mat',v_sim))

for iies = 1: length(Ies)
  iies
  for iiis = 1: length(Iis)
    for ei = 1 : length(EIs)
%       for igain = 1:length(Gains)
        %         igain
        load(sprintf('~/pmod/proc/ei/pmod_wc_ei_Ie%d_Ii%d_ei%d_v%d.mat',iies,iiis,ei,v))

        if round(Ies(iies)*100)/100 == -2.8 && round( Iis(iiis)*100)/100 == -3.4
          disp('save stuff')
          FFCC = out.FC;
        end
        
        % Time scales
        outp.lambda(:,iies,iiis,ei,igain)      = mean(out.lambda,2);
        outp.lambda_env(:,iies,iiis,ei,igain)  = mean(out.lambda_env,2);
        % DFA
        outp.dfa_sim(:,iies,iiis,ei,igain)     = squeeze(mean(out.dfa,1));
        outp.dfa_env_sim(:,iies,iiis,ei,igain) = squeeze(mean(out.dfa_env,1));
        
        outp.peakfreq(iies,iiis,ei,igain) = out.peakfreq;
        
        pars.dim = 2;
        
        outp.fc_sim_tmp = tp_match_aal(pars,out.FC,3);  
        outp.fc_sim_env_tmp = tp_match_aal(pars,out.FC_env,3);  
                
        [outp.r_rest_corr(iies,iiis,ei,igain), outp.p_rest_corr(iies,iiis,ei,igain)]=corr(outp.fc_sim_tmp(mask),fc_rest(mask));
        [outp.r_env_rest_corr(iies,iiis,ei,igain), outp.p_env_rest_corr(iies,iiis,ei,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_rest(mask));
        [outp.r_env_task_corr(iies,iiis,ei,igain), outp.p_env_task_corr(iies,iiis,ei,igain)]=corr(outp.fc_sim_env_tmp(mask),fc_task(mask));

%         [outp.r_rest_corr_avg(iies,iiis,iG,igain), outp.p_rest_corr_avg(iies,iiis,iG,igain)]=corr(nanmean(outp.fc_sim_tmp)',nanmean(fc_rest)'); 
%         outp.r_rest_corr_unc(iies,iiis,iG,igain) = dot(outp.fc_sim_tmp(mask),fc_rest(mask)) / sqrt(dot(outp.fc_sim_tmp(mask),outp.fc_sim_tmp(mask)) * dot(fc_rest(mask),fc_rest(mask)));
   
        pars.dim = 1;
        
        [outp.dfa_r_rest(iies,iiis,ei,igain), outp.dfa_p_rest(iies,iiis,ei,igain)] = corr(dfa_emp_rest(:),tp_match_aal(pars,repmat(outp.dfa_sim(:,iies,iiis,ei,igain),[1 90]),pars));
       	[outp.dfa_env_r_rest(iies,iiis,ei,igain), outp.dfa_env_p_rest(iies,iiis,ei,igain)] = corr(dfa_emp_rest(:),tp_match_aal(pars,repmat(outp.dfa_env_sim(:,iies,iiis,ei,igain),[1 90]),pars));

        [outp.lambda_r_rest(iies,iiis,ei,igain), outp.lambda_p_rest(iies,iiis,ei,igain)] = corr(lambda_emp_rest,tp_match_aal(pars,repmat(outp.lambda(:,iies,iiis,ei,igain),[1 90]),pars));
        [outp.lambda_r_task(iies,iiis,ei,igain), outp.lambda_p_task(iies,iiis,ei,igain)] = corr(lambda_emp_task,tp_match_aal(pars,repmat(outp.lambda(:,iies,iiis,ei,igain),[1 90]),pars));
        
        [outp.lambda_env_r_rest(iies,iiis,ei,igain), outp.lambda_env_p_rest(iies,iiis,ei,igain)] = corr(lambda_emp_rest,tp_match_aal(pars,repmat(outp.lambda_env(:,iies,iiis,ei,igain),[1 90]),pars));

        outp.dist_fc_rest (iies,iiis,ei,igain)  = mean(outp.fc_sim_tmp(mask))-mean(fc_rest(mask));
        outp.dist_fc_task (iies,iiis,ei,igain)  = mean(outp.fc_sim_tmp(mask))-mean(fc_task(mask));
        
        outp.fc_sim_mean(iies,iiis,ei,igain) = mean(outp.fc_sim_tmp(mask));
        outp.fc_sim_env_mean(iies,iiis,ei,igain) = mean(outp.fc_sim_env_tmp(mask));
        
        % KURAMOTO
        outp.kuramoto_mean (iies,iiis,ei,igain) = mean(out.KOPmean);
        outp.kuramoto_std (iies,iiis,ei,igain)  = mean(out.KOPsd);
        
        outp.kuramoto_mean_diff (iies,iiis,ei,igain) = mean(out.KOPmean) - mean(kura_mean(:,1,1,1));
        outp.kuramoto_std_diff (iies,iiis,ei,igain)  = mean(out.KOPsd)- mean(kura_std(:,1,1,1));
        
        fclose all;
        %
      end
    end
end

  save(sprintf('~/pmod/proc/pmod_wc_ei_all_v%d.mat',v_sim),'outp')
% % % 
else
  load(sprintf('~/pmod/proc/pmod_wc_ei_all_v%d.mat',v_sim))
end

error('!')

%% PLOT GAIN VS NO GAIN
% increase in gain vs. baseline
clear par
load(sprintf('~/pmod/proc/pmod_final_fitting_fits_v%d.mat',4))

idx = [find(Ies==par.rest(1)) find(Iis==par.rest(2))];
idx2 = [find(Ies==par.task(1)) find(Iis==par.task(2))];
[i,idx3(:,1)]=find(Ies==par.task_alt(:,1));
[i,idx3(:,2)]=find(Iis==par.task_alt(:,2));

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
 	scatter(ax{iax},idx(2),idx(1),20,'markerfacecolor','w','markeredgecolor','k')
  scatter(ax{iax},idx2(2),idx2(1),20,'markerfacecolor','y','markeredgecolor','k')
  scatter(ax{iax},idx3(:,2),idx3(:,1),20,'markerfacecolor','y','markeredgecolor','k')
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


% BAR PLOTS
diag_mask = diag(ones(6,6));


par1(:,iei) = squeeze(abs(outp.fc_sim_env_mean(idx(1),idx(2),iei)))-abs(squeeze(outp.fc_sim_env_mean(idx(1),idx(2),7)))
par2(:,iei) = squeeze(abs(outp.fc_sim_env_mean(idx2(1),idx2(2),iei)))-abs(squeeze(outp.fc_sim_env_mean(idx2(1),idx2(2),7)))
par3 = squeeze(abs(outp.fc_sim_env_mean(idx3(:,1),idx3(:,2),iei)))-abs(squeeze(outp.fc_sim_env_mean(idx3(:,1),idx3(:,2),7)))
all_par3(:,iei) = diag(par3);

figure;set(gcf,'color','w')
subplot(2,2,1)
bar([1 2 3 4 5 6 7 8],[par1(:,iei); par2(:,iei); all_par3(:,iei)])
tp_editplots; axis square
axis([0 9 -0.03 0.03])

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_ei_eivsbaseline_bar_ei%d_v%d.pdf',iei,v_sim))

%% 



