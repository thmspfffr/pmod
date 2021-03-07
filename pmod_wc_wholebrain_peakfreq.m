%% pmod_wc_wholebrain_peakfreq
% Stochastic simulation of 2*N WC nodes during "rest"
%--------------------------------------------------------------------------

clear
% restoredefaultpath;matlabrc

outdir = '~/pmod/proc/';

% 29-05-2018: fit E and I through resting state recordings. Then obtain
% changes in E and I due to task from recordings and keep those parameters
% fixed for the drug simulations. Vary excitability and gain for the drug
% recordings.
%--------------------------------------------------------------------------

%-------------------------------------------------------------------------
% VERSION 1: 20-10-2018
% %-------------------------------------------------------------------------
% v           = 1;
% Ies         = -4:0.05:-1;
% Iis         = -5:0.05:-2;
% Gg          = 0:0.05:2;
% Gains       = [0 0.45];
% nTrials     = 1;
% tmax        = 6500;  % in units of tauE
%-------------------------------------------------------------------------
% VERSION 1: 20-10-2018
% %-------------------------------------------------------------------------
% v           = 2;
% Ies         = -4:0.025:-1;
% Iis         = -5:0.025:-2;
% Gg          = 1.4;
% Gains       = [0:0.05:0.6 -0.2:0.05:-0.05 0.65:0.05:1];
% nTrials     = 1;
% tmax        = 6500;  % in units of tauE
%-------------------------------------------------------------------------
% VERSION 3: M
% %-------------------------------------------------------------------------
v           = 3;
Ies         = -4:0.025:-1;
Iis         = -5:0.025:-2;
Gg          = [1.2:-0.01:1.10];
Gains       = [-0.1:0.02:0.12];
nTrials     = 5;
tmax        = 6500;  % in units of tauE
%-------------------------------------------------------------------------

dt          = 0.01;

% EXCLUDE CERTAIN REGIONS - BCN ordering
k = 1 : 90;
exclude_bcn = [11 15 21 36 37 38 39 52 53 54 55 70 76 80];
include_bcn = find(~ismember(k,exclude_bcn));
% ------------------

load ~/sc90.mat %Bea SC
C = SC/max(SC(SC>0));
C = C(include_bcn,include_bcn);
N = size(C,1);

%--------------------------------------------------------------------------
% PARAMETER DEFINITIONS
%--------------------------------------------------------------------------
% Connectivity:
wII=4;
wIE=16;
wEI=12;
wEE=12;
% noise level
sigma = 0.0005;

tauE = 1;
tauI = 2;
tau = zeros(2*N,1);
tau(1:N) = tauE;
tau(N+1:2*N) = tauI;


tspan=0:dt:tmax;
L = length(tspan);
clear tspan

ds = 10;
Tds = length(0:ds*dt:tmax)-1;
tauEsec = 0.009; % in seconds
resol = ds*dt*tauEsec;

isub = find( triu(ones(N)) - eye(N) );

% ----------------
% WAVELETS
% ----------------
% [wavelet,f]=tp_mkwavelet(11.3137,0.5,(1/resol));
% delta_time = 6/pi./(f(2)-f(1));
% delta_time = round(delta_time*1000)/1000;
% t_shift    = delta_time;
% n_win      = round(delta_time*(1/resol));
% n_shift    = round(t_shift*(1/resol));
% nseg=floor((L/10-n_win)/n_shift+1);
% ----------------

% LOAD REGIONS THAT ARE SIGNIFICANTLY (P<0.05) ALTERED BY TASK (mostly
% decreased), in the relevant freq range
vvv = 33; % AAL
load(sprintf('~/pupmod/proc/src/pupmod_src_pow_taskmodulationindex_v%d.mat',vvv))


% load(sprintf('~/pmod/proc/pmod_final_fitting_indivfits_rest_v%d.mat',3))
%
% clear idx

% task_exc = 8;
% task_inh = 15;

% task_exc = 10;
% task_inh = 17;
%
% [Mu,ia,ic]=unique([idx_rest.inh' idx_rest.exc'],'rows','stable')
% c=accumarray(ic,1);
% cols = cbrewer('seq', 'YlOrRd', max(c)+3);
% cols = cols(4:end,:);
%
% % find idx of delta_gain = 0 and coupling (Gg) = 1.15
% baseline_gain=find(Gains==0);
% baseline_coupling = find(Gg==1.15);
%
% SUBJ = 1:28;
% clear fc_mod_rest fc_mod_task

load(sprintf('~/pmod/proc/pmod_final_fitting_indivfits_rest_v%d.mat',v))

% taken from final_fitting
task_exc = 10;
task_inh = 17;

task_Es = 0;
task_Is = 0;

%%

for isubj = 1:size(idx_rest.exc,2)
  
  
  
  if ~exist(sprintf('~/pmod/proc/numerical/peakfreq/v%d/',v))
    mkdir(sprintf(['~/pmod/proc/numerical/peakfreq/' 'v%d'],v))
  end
  
  % save configuration, with all above parameters
  if ~exist(sprintf([outdir 'pmod_wc_wholebrain_peakfreq_parameters_v%d.mat'],v))
    save(sprintf([outdir 'pmod_wc_wholebrain_peakfreq_parameters_v%d.mat'],v))
  end
  
  outdir = sprintf(['~/pmod/proc/numerical/peakfreq/v%d/'],v);
  
  fn = sprintf('pmod_wc_wholebrain_peakfreq_isubj%d_v%d',isubj,v);
  if tp_parallel(fn,outdir,1,0)
    continue
  end
%   
  for icond = 1 : 2
    
    if icond == 1
      iies = idx_rest.exc(isubj);
      iiis = idx_rest.inh(isubj);
    else
      iies = idx_rest.exc(isubj)+task_exc;
      iiis = idx_rest.inh(isubj)+task_inh;  
    end
    
    for igain = 1:length(Gains)
      for iG = 1:length(Gg)
        itask_E = 1 ;
        itask_I = 1 ;

        fprintf('Subject%d, Gain%d, Coupl%d...\n',isubj,igain,iG)
        
        g = Gg(iG);
        W = [wEE*eye(N)+g*C -wEI*eye(N); wIE*eye(N) -wII*eye(N)];
        
        % transfer function:
        gE = 1 + Gains(igain);
        gI = 1 + Gains(igain);
        aE  = 1/gE;
        aI  = 1/gI;
        Fe  = @(x) 1./(1 + exp(-x/aE) );
        Fi  = @(x) 1./(1 + exp(-x/aI) );
        F   = @(x) [feval(Fe,x(1:N));feval(Fi,x(N+1:2*N))];
        
        % Working point:
        Io          = zeros(2*N,1);
        Io(1:N)     = Ies(iies);
        Io(1:N)     = Io(1:N)+task_idx.*task_Es(itask_E); % add task effect to E
        Io(N+1:2*N) = Iis(iiis);
        Io(N+1:2*N) = Io(N+1:2*N)+task_idx.*task_Is(itask_I); % add task effect to I
        
        T           = Tds*resol; %% define time of interval
        freqs       = (0:Tds/2)/T; %% find the corresponding frequency in Hz
        freq100     = freqs(freqs<100 & freqs>1);
        pp          = 1:10:length(freq100);
        PSD         = zeros(length(pp),N,nTrials);
        frequencies = freq100(1:10:end)';
        
        % RUN SIMULATION
        % ---------------------
        
        for tr=1:nTrials
          r   = 0.001*rand(2*N,1);
          R   = zeros(Tds,N);
          Ri  = zeros(Tds,N);
          tt  = 0;
          %         transient:
          for t = 1:20000
            u = W*r + Io;
            K = feval(F,u);
            r = r + dt*(-r + K)./tau + sqrt(dt)*sigma*randn(2*N,1);
          end
          %         simulation
          for t = 1:L
            %           100*t/L
            u = W*r + Io;
            K = feval(F,u);
            r = r + dt*(-r+ K)./tau + sqrt(dt)*sigma*randn(2*N,1);
            if mod(t,ds)==0
              tt=tt+1;
              R(tt,:)  = r(1:N);
              Ri(tt,:)  = r(N+1:end);
            end
          end
          
          rE = R;
          rI = Ri;
          clear R Ri
          
          % FC matrix
          % ---------------------
          rc       	= corrcoef(rE);
          covariance = cov(rE);
%           covariance(find(eye(size(covariance,1))))=1;
          out.hc(igain,iG,icond,tr) = diff_entropy(covariance);
          fc      	= rc(isub);
          out.FR(:,igain,iG,icond,tr) = mean(rE);
          out.fc_FR(:,:,igain,iG,icond,tr) = rc;
          
          for i=1:N
            
            % COMPUTE POWER SPECTRUM
            % ---------------------------
            f = rE(:,i) - mean(rE(:,i));
            xdft = fft(f);
            xdft = xdft(1:floor(Tds/2)+1);
            pw = (1/(Tds/2)) * abs(xdft).^2;
            psd = pw(freqs<100 & freqs>1);
            f = freqs(freqs<100 & freqs>1);
            fnew = f(1:10:end);
            psd  = psd(1:10:end);
            PSD(:,i,tr) = psd';
            f = fnew;
            
          end
          % ----------------------------------
          % EXTRACT PEAK FREQ
          % ---------------------------
          [~,peak_idx] = max(smooth(mean(PSD(f>4,:),2),20));
          out.peakfreq(igain,iG,icond) = f(peak_idx+find(f<4,1,'last'));
          clear PSD rc fc_env
          
          z                                 = rE + 1i*rI;
          ku                                = sum(z,2)/N;
          out.KOP                           = abs(ku);
          out.KOPsd(igain,iG,icond,tr,1)    = std(out.KOP);
          out.KOPmean(igain,iG,icond,tr,1)  = mean(out.KOP);
          
        end
        
      end
    end
  end
  
  save(sprintf([outdir '/%s.mat'],fn),'out')

  % make sure file is saved
  while 1
    try
      load(sprintf([outdir '%s.mat'],fn))
      break
    catch me
      save(sprintf([outdir '%s.mat'],fn),'out')
    end
  end
  tp_parallel(fn,outdir,0,0)
end

exit
%% PLOT PEAK FREQ

igain_atx = 11; iG_atx = 6
igain_dpz = 8; iG_dpz = 10

v = 3
mask = logical(tril(ones(76,76),-1));

for isubj = 1 : 28
  
    load(sprintf('~/pmod/proc/numerical/peakfreq/v%d/pmod_wc_wholebrain_peakfreq_isubj%d_v%d.mat',v,isubj,v))
    tmp1 = squeeze(nanmean(out.fc_FR(:,:,6,6,1,:),6));
    fc_model(:,isubj,1)=tmp1(mask);
    tmp2 = squeeze(nanmean(out.fc_FR(:,:,6,6,2,:),6));
    fc_model(:,isubj,2)=tmp2(mask);
%   pf(:,:,:,isubj) = out.peakfreq;
%   KOPsd(isubj,:,:,:)=out.KOPsd;
%   KOPmean(isubj,:,:,:)=out.KOPmean;
%   hc(:,:,:,isubj) = out.hc;
end

figure_w

m_rest = nanmean(nanmean(fc(:,:,1,1),1),2);
s_rest = nanstd(nanmean(fc(:,:,1,1),1),[],2)/sqrt(28);
m_task = nanmean(nanmean(fc(:,:,1,2),1),2);
s_rest = nanstd(nanmean(fc(:,:,1,2),1),[],2)/sqrt(28);
prct = 100*(m_task-m_rest)/m_rest;

m_m_rest = nanmean(fc_model(:,1),1);
s_m_rest = nanstd(fc_model(:,1),[],1)/sqrt(28);
m_m_task = nanmean(fc_model(:,2),1);
s_m_rest = nanstd(fc_model(:,2),[],1)/sqrt(28);
prct_m = 100*(m_m_task-m_m_rest)/m_m_rest;

[h,~,~,s]=ttest(fc(:,:,1,2),fc(:,:,1,1),'dim',2);

altered_corr_neg = 100*sum(h>0 & s.tstat < 0)/size(fc,1);
altered_corr_pos = 100*sum(h>0 & s.tstat > 0)/size(fc,1);

[h,~,~,s]=ttest(fc_model(:,:,2),fc_model(:,:,1),'dim',2);

altered_corr_m_neg = 100*sum(h>0 & s.tstat < 0)/size(fc,1);
altered_corr_m_pos = 100*sum(h>0 & s.tstat > 0)/size(fc,1);

subplot(2,2,1)
bar([1 2],[m_rest m_task])
axis square; tp_editplots
axis([0.5 2.5 0 0.06])

subplot(2,2,2)
bar([1 2],[m_m_rest m_m_task])
axis square; tp_editplots
axis([0.5 2.5 0 0.04])

% 
% figure_w;
% subplot(2,2,1);
% imagesc(hc(:,:,1,1),[-5 5]); colormap(cmap); axis square
% set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6); xlabel('Change in gain')
% set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6); ylabel('Change in coupling')
% tp_editplots
% subplot(2,2,2);
% imagesc(hc(:,:,1,2),[-5 5]); colormap(cmap); axis square
% set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6); xlabel('Change in gain')
% set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6); ylabel('Change in coupling')
% tp_editplots
% subplot(2,2,3);
% imagesc(hc(:,:,1,3),[-5 5]); colormap(cmap); axis square
% set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6); xlabel('Change in gain')
% set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6); ylabel('Change in coupling')
% tp_editplots
% subplot(2,2,4);
% imagesc(hc(:,:,1,4),[-5 5]); colormap(cmap); axis square
% set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6); xlabel('Change in gain')
% set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6); ylabel('Change in coupling')
% tp_editplots
% 
% print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_entropy_v%d.pdf',v))
%%


figure_w
k=subplot(2,2,1); hold on
imagesc(squeeze(nanmean(KOPmean(:,:,:,1),1))-squeeze(nanmean(KOPmean(:,6,6,1),1)),[-0.05 0.05])
plot(6,6,'o','markeredgecolor','k','markerfacecolor','w','markersize',5)
colormap(k,cmap); axis square tight; tp_editplots; set(gca,'ydir','reverse')
plot(igain_atx,iG_atx,'o','markeredgecolor','k','markerfacecolor','r','markersize',5)
plot(igain_dpz,iG_dpz,'o','markeredgecolor','k','markerfacecolor','b','markersize',5)
set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6); xlabel('Change in gain')
set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6); ylabel('Change in coupling')

k=subplot(2,2,2); hold on
imagesc(squeeze(nanmean(KOPmean(:,:,:,2),1))-squeeze(nanmean(KOPmean(:,6,6,2),1)),[-0.05 0.05])
plot(6,6,'o','markeredgecolor','k','markerfacecolor','w','markersize',5)
colormap(k,cmap); axis square tight; tp_editplots; set(gca,'ydir','reverse')
plot(igain_atx,iG_atx,'o','markeredgecolor','k','markerfacecolor','r','markersize',5)
plot(igain_dpz,iG_dpz,'o','markeredgecolor','k','markerfacecolor','b','markersize',5)
set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6); xlabel('Change in gain')
set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6); ylabel('Change in coupling')

k=subplot(2,2,3); hold on
imagesc(squeeze(nanmean(KOPsd(:,:,:,1),1))-squeeze(nanmean(KOPsd(:,6,6,1),1)),[-0.0002 0.0002])
plot(6,6,'o','markeredgecolor','k','markerfacecolor','w','markersize',5)
colormap(k,cmap); axis square tight; tp_editplots; set(gca,'ydir','reverse')
plot(igain_atx,iG_atx,'o','markeredgecolor','k','markerfacecolor','r','markersize',5)
plot(igain_dpz,iG_dpz,'o','markeredgecolor','k','markerfacecolor','b','markersize',5)
set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6); xlabel('Change in gain')
set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6); ylabel('Change in coupling')

k=subplot(2,2,4); hold on
imagesc(squeeze(nanmean(KOPsd(:,:,:,2),1))-squeeze(nanmean(KOPsd(:,6,6,2),1)),[-0.0002 0.0002])
plot(6,6,'o','markeredgecolor','k','markerfacecolor','w','markersize',5)
colormap(k,cmap); axis square tight; tp_editplots; set(gca,'ydir','reverse')
plot(igain_atx,iG_atx,'o','markeredgecolor','k','markerfacecolor','r','markersize',5)
plot(igain_dpz,iG_dpz,'o','markeredgecolor','k','markerfacecolor','b','markersize',5)
set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6); xlabel('Change in gain')
set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6); ylabel('Change in coupling')

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_kuramoto_v%d.pdf',v))


figure_w

k=subplot(2,2,1); hold on
imagesc(nanmean(pf(:,:,1,:),4),[7 13]); colormap(cmap)
plot(6,6,'o','markeredgecolor','k','markerfacecolor','w','markersize',5)
plot(igain_atx,iG_atx,'o','markeredgecolor','k','markerfacecolor','r','markersize',5)
plot(igain_dpz,iG_dpz,'o','markeredgecolor','k','markerfacecolor','b','markersize',5)
colormap(k,plasma); axis square tight; tp_editplots; set(gca,'ydir','reverse')
set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6); xlabel('Change in gain')
set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6); ylabel('Change in coupling')

l=subplot(2,2,2); hold on
imagesc(nanmean(pf(:,:,2,:),4),[7 13]); colormap(cmap)
plot(6,6,'o','markeredgecolor','k','markerfacecolor','w','markersize',5)
plot(igain_atx,iG_atx,'o','markeredgecolor','k','markerfacecolor','r','markersize',5)
plot(igain_dpz,iG_dpz,'o','markeredgecolor','k','markerfacecolor','b','markersize',5)
axis square tight; tp_editplots; set(gca,'ydir','reverse')
set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6); xlabel('Change in gain')
set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6); ylabel('Change in coupling')

subplot(2,2,3); hold on
imagesc(nanmean(pf(:,:,1,:),4)-nanmean(pf(6,6,1,:),4),[-2 2]); colormap(cmap)
plot(6,6,'o','markeredgecolor','k','markerfacecolor','w','markersize',5)
plot(igain_atx,iG_atx,'o','markeredgecolor','k','markerfacecolor','r','markersize',5)
plot(igain_dpz,iG_dpz,'o','markeredgecolor','k','markerfacecolor','b','markersize',5)
colormap(cmap); axis square tight; tp_editplots; set(gca,'ydir','reverse')
set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6); xlabel('Change in gain')
set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6); ylabel('Change in coupling')


subplot(2,2,4); hold on
imagesc(nanmean(pf(:,:,2,:),4)-nanmean(pf(6,6,2,:),4),[-2 2]); colormap(cmap)
plot(6,6,'o','markeredgecolor','k','markerfacecolor','w','markersize',5)
plot(igain_atx,iG_atx,'o','markeredgecolor','k','markerfacecolor','r','markersize',5)
plot(igain_dpz,iG_dpz,'o','markeredgecolor','k','markerfacecolor','b','markersize',5)
colormap(cmap); axis square tight; tp_editplots; set(gca,'ydir','reverse')
set(gca,'xtick',[1 6 11],'xticklabel',[-0.1 0 0.1],'fontsize',6); xlabel('Change in gain')
set(gca,'ytick',[1 6 11],'yticklabel',[0.05 0 -0.05],'fontsize',6); ylabel('Change in coupling')

colormap(l,plasma); 
colormap(k,plasma); 

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_peakfreq_v%d.pdf',v))

%% EMPIRICAL KURA
SUBJLIST        = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];
ord    = pconn_randomization;
kura_std = []; kura_mean = [];

for isubj = 1:length(SUBJLIST)
  for m = 1 : 3
    im = find(ord(SUBJLIST(isubj),:)==m);
    load(sprintf('/home/tpfeffer/pupmod/proc/conn/pupmod_src_kuramoto_s%d_m%d_v%d.mat',SUBJLIST(isubj),im,v))
    kura_std(isubj,m,:) = nanmean(kuramoto.Rsd,2);
    kura_mean(isubj,m,:) = nanmean(kuramoto.Rmean,2);
    try
    load(sprintf('/home/tpfeffer/pupmod/proc/conn/pupmod_task_src_kuramoto_s%d_m%d_v%d.mat',SUBJLIST(isubj),im,v))
    kura_std_task(isubj,m,:) = nanmean(kuramoto.Rsd,2);
    kura_mean_task(isubj,m,:) = nanmean(kuramoto.Rmean,2);
    
    catch me
      kura_std_task(isubj,m,:) = nan(17,1);
    kura_mean_task(isubj,m,:) = nan(17,1);
    end
      
  end
end

cols = [0.7 0.7 0.7; 1 0 0; 0 0 1]


figure_w
subplot(2,2,1); hold on
plot(log10(freqoi),squeeze(nanmean(kura_mean))')
set(gca,'xtick',log10(freqoi([1 5 9 13 17])),'xticklabel',freqoi([1 5 9 13 17]),'fontsize',6)
tp_editplots; xlabel('Frequency [Hz]'); ylabel('Synchrony (R)')

subplot(2,2,2); hold on
plot(log10(freqoi),squeeze(nanmean(kura_std))')
set(gca,'xtick',log10(freqoi([1 5 9 13 17])),'xticklabel',freqoi([1 5 9 13 17]),'fontsize',6)
tp_editplots; xlabel('Frequency [Hz]'); ylabel('Metastability (R_sd)')

subplot(2,2,3); hold on
plot(log10(freqoi),squeeze(nanmean(kura_mean_task))')
set(gca,'xtick',log10(freqoi([1 5 9 13 17])),'xticklabel',freqoi([1 5 9 13 17]),'fontsize',6)
tp_editplots; xlabel('Frequency [Hz]'); ylabel('Synchrony (R)')

subplot(2,2,4); hold on
plot(log10(freqoi),squeeze(nanmean(kura_std_task ))')
set(gca,'xtick',log10(freqoi([1 5 9 13 17])),'xticklabel',freqoi([1 5 9 13 17]),'fontsize',6)
tp_editplots; xlabel('Frequency [Hz]'); ylabel('Metastability (R_sd)')

colororder(cols);

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_kuramoto_meg_v%d.pdf',v))
