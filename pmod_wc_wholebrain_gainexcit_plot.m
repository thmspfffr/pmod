%% FITTING
% pmod_wc_wholebrain_plot.m
%-------------------------------------------------------------------------
% Version 1 - time is shorter than usual
%-------------------------------------------------------------------------
v_sim           = 2;
Ies         = [-2.8 -1.8];
Iis         = [-3.5 -2.4];
Gg          = 0.62;
Gains       = -0.25:0.025:0.25;
Excit       = 0.8:0.025:1.2;
nTrials     = 3;
tmax        = 65000; % in units of tauE
% wins        = [3 50];
%-------------------------------------------------------------------------

% connectivity, AAL
v_conn =  1;
% simulations
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
      igain
      for iexc = 1 : length(Excit)
        %         igain
        load(sprintf('~/pmod/proc/pmod_wc_wholebrain_gainexcit_ei%d_gain%d_excit%d_v%d.mat',ei,igain,iexc,v_sim))
       
        
        % Time scales
        outp.lambda(:,ei,igain,iexc)      = mean(out.lambda,2);
        outp.lambda_env(:,ei,igain,iexc)  = mean(out.lambda_env,2);
        % DFA
%         outp.dfa_sim(:,iies,iiis,iG,igain)     = squeeze(mean(out.dfa,1));
%         outp.dfa_env_sim(:,iies,iiis,iG,igain) = squeeze(mean(out.dfa_env,1));
        
        [~,peak_idx]=max(smooth(mean(out.PSD(out.f>3,:),2),20));
        outp.peakfreq(ei,igain,iexc) = out.f(peak_idx+find(out.f<3,1,'last'));
        
        pars.dim = 2;
        outp.fc_sim_tmp = tp_match_aal(pars,mean(out.FC,3));
        outp.fc_sim_var(:,ei,igain,iexc) = std(nanmean(outp.fc_sim_tmp)./max(nanmean(outp.fc_sim_tmp)));
        
        fc_sim_env = tp_match_aal(pars,mean(out.FC_env,3));
        outp.fc_env_sim_mean(ei,igain,iexc) = mean(fc_sim_env(mask));
       
        [outp.r_rest_corr(ei,igain,iexc), outp.p_rest_corr(ei,igain,iexc)]=corr(outp.fc_sim_tmp(mask),fc_rest(mask));
        [outp.r_rest_corr_avg(ei,igain,iexc), outp.p_rest_corr_avg(ei,igain,iexc)]=corr(nanmean(outp.fc_sim_tmp)',nanmean(fc_rest)');
        [outp.r_env_rest_corr(ei,igain,iexc), outp.p_env_rest_corr(ei,igain,iexc)]=corr(fc_sim_env(mask),fc_rest(mask));

        outp.r_rest_corr_unc(ei,igain,iexc) = dot(outp.fc_sim_tmp(mask),fc_rest(mask)) / sqrt(dot(outp.fc_sim_tmp(mask),outp.fc_sim_tmp(mask)) * dot(fc_rest(mask),fc_rest(mask)));
        
        pars.dim = 1;

        %
        [outp.lambda_r_rest(ei,igain,iexc), outp.lambda_p_rest(ei,igain,iexc)] = corr(lambda_emp_rest,tp_match_aal(pars,repmat(outp.lambda(:,ei,igain,iexc),[1 90]),pars));
        [outp.lambda_r_task(ei,igain,iexc), outp.lambda_p_task(ei,igain,iexc)] = corr(lambda_emp_task,tp_match_aal(pars,repmat(outp.lambda(:,ei,igain,iexc),[1 90]),pars));
        [outp.lambda_env_r_rest(ei,igain,iexc), outp.lambda_env_p_rest(ei,igain,iexc)] = corr(lambda_emp_rest,tp_match_aal(pars,repmat(outp.lambda_env(:,ei,igain,iexc),[1 90]),pars));
       
        outp.dist_fc_rest (ei,igain,iexc)  = mean(outp.fc_sim_tmp(mask))-mean(fc_rest(mask));
        outp.dist_fc_task (ei,igain,iexc)  = mean(outp.fc_sim_tmp(mask))-mean(fc_task(mask));
        
        outp.fc_sim_mean(ei,igain,iexc) = mean(outp.fc_sim_tmp(mask));
        outp.fc_sim_all(:,:,ei,igain,iexc) = outp.fc_sim_tmp;
        %         outp.fc_sim_env(iies,iiis,iG,igain) = mean(outp.fc_sim_env_tmp(mask));
        %
        outp.Ies(iies) = out.Ie;
        outp.Iis(iiis) = out.Ii;
        
        % KURAMOTO
        outp.kuramoto_mean (ei,igain,iexc) = mean(out.KOPmean);
        outp.kuramoto_std (ei,igain,iexc)  = mean(out.KOPsd);
        
        fclose all;
        %
      end
    end
  end
  
  
  % -----------------------------------
  % COMPUTE DEGREE IN MODEL
  % -----------------------------------
  
%   vec = 1 : 90;
%   
%   for ei = 1 : 2
%     ei
%     for igain = 1 : length(Gains)
%       for iexc = 1 : length(Excit)
%         
%         for ij = 1 : 90
%           for jj = 1 : 90
%             
%             jjj = vec(vec~=ij);
%             fc_tmp = outp.fc_sim_all(ij,jj,ei,igain,iexc);
%             
%             x_ref1 = mean(outp.fc_sim_all(ij,jjj,ei,igain,iexc),2);
%             
%             th(ij,jj,ei,igain,iexc) = fc_tmp>x_ref1;
%           end
%         end
%         th_all(:,:,ei,igain,iexc) = th(:,:,ei,igain,iexc)+fliplr(rot90(th(:,:,ei,igain,iexc)>0,-1));
%       end
%     end
%     
%     for igain = 1 : length(Gains)
%       for ii = 1 : length(Iis)
%         for ie = 1 : length(Ies)
%           outp.degree(ei,igain,iexc) = nanmean(nanmean(th_all(:,:,ei,igain,iexc)));
%         end
%       end
%     end
%   end
  
  save(sprintf('~/pmod/proc/pmod_wholebrain_gainexcit_all_v%d.mat',v_sim),'outp')
  
else
  load(sprintf('~/pmod/proc/pmod_wholebrain_gainexcit_all_v%d.mat',v_sim))
end

error('!')

%% PLOT BASIC PARAMETERS
% Lambda
% Peak frequency
% Correlation model with MEG: Lambda and FC
% ------------------------------------------
clear par
ei = 1
% idx = [find(round(Ies*100)/100==-2.8) find(round(Iis*100)/100==-3.5000)];
% idx2 = [find(round(Ies*100)/100==-1.8) find(round(Iis*100)/100==-2.4000)];

% igain = 3;
% G = 1;

figure; set(gcf,'color','w')

% plot lambda
ax{1} = subplot(2,2,1); hold on
par = squeeze(outp.fc_sim_mean(ei,:,:));
% par(osc1(:,:,:,3)>0.5)=nan;
imagescnan(par,[0 0.02])
title('FC_{FR}');

% plot correlation FC model / MEG
ax{2} = subplot(2,2,2); hold on
par = squeeze(outp.fc_env_sim_mean(ei,:,:));
% par(osc1(:,:,:,3)>0.5)=nan;
imagescnan(par,[0 0.02])
title('FC_{Env}');

% plot correlation lambda model / MEG
ax{3} = subplot(2,2,3); hold on
par = squeeze(outp.r_rest_corr(ei,:,:));
% par(osc1(:,:,:,3)>0.5)=nan;
imagescnan(par,[0 0.4])
title('r(FC_{FR})');

% plot peak freq model
ax{4} = subplot(2,2,4); hold on
par = squeeze(outp.r_env_rest_corr(ei,:,:));
% par(osc1(:,:,:,3)>0.5)=nan;
imagescnan(par,[0 0.4])
title('r(FC_{Env})');

for iax = 1 : length(ax)
  
  if iax == 1
    ylabel(ax{iax},'Gains');
    set(ax{iax},'YTick',1:5:length(Gains ),'YTickLabels',Gains (1:5:end))
    set(ax{iax},'XTick',1:4:length(Excit),'XTickLabels',Excit (1:4:end))
  elseif iax == 2
    set(ax{iax},'YTick',1:5:length(Gains ),'YTickLabels',Gains (1:5:end))
    set(ax{iax},'XTick',1:4:length(Excit),'XTickLabels',Excit (1:4:end))
  elseif iax == 3
    xlabel(ax{iax},'Excitability')
    ylabel(ax{iax},'Gains');
    set(ax{iax},'YTick',1:5:length(Gains ),'YTickLabels',Gains (1:5:end))
    set(ax{iax},'XTick',1:4:length(Excit),'XTickLabels',Excit (1:4:end))
  elseif iax ==4
    xlabel(ax{iax},'Excitability')
    set(ax{iax},'YTick',1:5:length(Gains ),'YTickLabels',Gains (1:5:end))
    set(ax{iax},'XTick',1:4:length(Excit),'XTickLabels',Excit (1:4:end))
  end
  tp_editplots(ax{iax})
  colormap(plasma)
  c = colorbar(ax{iax}); axis(ax{iax},'tight')
  c.Ticks = [min(c.Ticks) max(c.Ticks)];
 
end

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_lambda_excitgain_ei%d_v%d.pdf',ei,v_sim))

%% PLOT GAIN VS NO GAIN
% increase in gain vs. baseline
clear par

ei = 1;

cmap = cbrewer('div', 'RdBu', 100,'pchip');
cmap = cmap(end:-1:1,:);

figure; set(gcf,'color','w')

% plot correlation FC model / MEG
ax{1} = subplot(2,2,1); hold on
par = squeeze(abs(outp.fc_sim_mean(ei,:,:)))-abs(squeeze(outp.fc_sim_mean(ei,11,9)));
% par(osc1>0.5)=nan;
imagescnan(par,[-0.01 0.01])
title('Contrast: FC_{FR}');

% plot correlation FC model / MEG
ax{2} = subplot(2,2,2); hold on
par = squeeze(abs(outp.fc_env_sim_mean(ei,:,:)))-abs(squeeze(outp.fc_env_sim_mean(ei,11,9)));
% par(osc1>0.5)=nan;
imagescnan(par,[-0.002 0.002])
title('Contrast: FC_{Env}');


% plot correlation lambda model / MEG
ax{3} = subplot(2,2,3); hold on
par = squeeze(mean(1./outp.lambda(:,ei,:,:)))-squeeze(mean(1./outp.lambda(:,ei,11,9)));
% par(osc1>0.5)=narn;
imagescnan(par,[-10 10])
title('Contrast: \lambda_{FR}');

% plot correlation lambda model / MEG
ax{4} = subplot(2,2,4); hold on
par = squeeze(mean(1./outp.lambda_env(:,ei,:,:)))-squeeze(mean(1./outp.lambda_env(:,ei,11,9)));
% par(osc1>0.5)=nan;
imagescnan(par,[-5 5])
title('Contrast: \lambda_{Env}');


for iax = 1 : length(ax)
  if iax == 1
    ylabel(ax{iax},'Gains');
    set(ax{iax},'YTick',1:5:length(Gains ),'YTickLabels',Gains (1:5:end))
    set(ax{iax},'XTick',1:4:length(Excit),'XTickLabels',Excit (1:4:end))
  elseif iax == 2
    set(ax{iax},'YTick',1:5:length(Gains ),'YTickLabels',Gains (1:5:end))
    set(ax{iax},'XTick',1:4:length(Excit),'XTickLabels',Excit (1:4:end))
  elseif iax == 3
    xlabel(ax{iax},'Excitability')
    ylabel(ax{iax},'Gains');
    set(ax{iax},'YTick',1:5:length(Gains ),'YTickLabels',Gains (1:5:end))
    set(ax{iax},'XTick',1:4:length(Excit),'XTickLabels',Excit (1:4:end))
  elseif iax ==4
    xlabel(ax{iax},'Excitability')
    set(ax{iax},'YTick',1:5:length(Gains ),'YTickLabels',Gains (1:5:end))
    set(ax{iax},'XTick',1:4:length(Excit),'XTickLabels',Excit (1:4:end))
  end
  tp_editplots(ax{iax})
  colormap(cmap)
  c = colorbar(ax{iax}); axis(ax{iax},'tight')
  c.Ticks = [min(c.Ticks) max(c.Ticks)];
 
end

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_wholebrain_gainvsbaseline_excitgain_ei%d_v%d.pdf',ei,v_sim))

