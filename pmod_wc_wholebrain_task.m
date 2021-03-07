%% pmod_wc_wholebrain_task
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
v           = 4;
Ies         = -4:0.025:-1;
Iis         = -5:0.025:-2;
Gg          = [1.2:-0.01:1.10];
Gains       = [-0.1:0.02:0.12];
nTrials     = 1;
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

load(sprintf('~/pmod/proc/pmod_final_fitting_indivfits_rest_v%d.mat',3))

task_Es = [-0.5:0.025:1];
task_Is = [-0.5:0.025:1.5];
% --------------------
% below are the task input levels from pmod_final_fitting.m
% --------------------
task_exc = 10;
task_inh = 17;
% --------------------
%%
for isubj = 1:size(idx_rest.exc,2)
  
  iies = idx_rest.exc(isubj);
  iiis = idx_rest.inh(isubj);
  
  if ~exist(sprintf('~/pmod/proc/numerical/task/v%d/',v))
    mkdir(sprintf(['~/pmod/proc/numerical/task/' 'v%d'],v))
  end
  
  % save configuration, with all above parameters
  if ~exist(sprintf([outdir 'pmod_wc_wholebrain_final_parameters_v%d.mat'],v))
    save(sprintf([outdir 'pmod_wc_wholebrain_final_parameters_v%d.mat'],v))
  end
  
  outdir = sprintf(['~/pmod/proc/numerical/task/v%d/'],v);
  
  fn = sprintf('pmod_wc_wholebrain_task_isubj%d_v%d',isubj,v);
  if tp_parallel(fn,outdir,1,0)
    continue
  end
  
  out.fc_FR = zeros(76,76,length(task_Es),length(task_Is));
  tic
  % REST (1) OR TASK (2)
  for icond = 1 : 2
  for igain = 1:length(Gains) % find(Gains==0)%
    for iG = 1:length(Gg) % find(Gg==1.15)

      toc
      fprintf('Subject%d, Gain%d, Coupl%d...\n',isubj,igain,iG)
      
      tic
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
      if icond == 1
        Io       = zeros(2*N,1);
        Io(1:N)  = Ies(iies);
        Io(N+1:2*N)= Iis(iiis);
      else
        Io          = zeros(2*N,1);
        Io(find(task_idx))  = Ies(iies+task_exc);
        Io(find(~task_idx)) = Ies(iies);

        Io(find(task_idx)+N) = Iis(iiis+task_inh);
        Io(find(~task_idx)+N)= Iis(iiis);
      end

      
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
        clear R Ri rI
        
        % FC matrix
        % ---------------------
        rc       	= corrcoef(rE);
        fc      	= rc(isub);
        
        out.fc_FR(:,:,igain,iG,icond) = rc;
        
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
        
        %           out.alphapow(:,tr) = squeeze(mean(PSD(frequencies>=flp&frequencies<=fhi,:,:),1));
        
        % ----------------------------------
        % WAVELET ANALYSIS
        % ----------------------------------
        %           for j=1:nseg
        %             dloc2=rE((j-1)*n_shift+1:(j-1)*n_shift+n_win,:)';
        %             dataf(j,:)=abs(dloc2*wavelet).^2;
        %           end
        %           out.rc_wl = corr(dataf);
        %           out.rc_wl_cov = cov(dataf);
        % ----------------------------------
        % EXTRACT PEAK FREQ
        % ---------------------------
        [~,peak_idx]=max(smooth(mean(PSD(f>4,:),2),20));
        out.peakfreq(igain,iG,icond) = f(peak_idx+find(f<4,1,'last'));
        clear PSD rc fc_env
        
      end
      
      
    end
    % end
  end
  end
  save(sprintf([outdir '/%s.mat'],fn),'out')
  %
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
error('!')
%%
v_sim = 4;

mask = logical(tril(ones(76,76),-1));
fc_all = zeros((76*76-76)/2, 12, 11,28);
for isubj = 1 : 28
  isubj
%   load(sprintf('~/pmod/proc/detosc/task/v%d/pmod_wc_wholebrain_detosc_task_isubj%d_v%d.mat',v_sim,isubj,v_sim))
%   osc(isubj,:,:) = squeeze(nanmean(out.osc1,1)>0.5);
%   osc(isubj,21,21)
  load(sprintf('~/pmod/proc/numerical/task/v%d/pmod_wc_wholebrain_task_isubj%d_v%d.mat',v_sim,isubj,v_sim))
  
  tmp_fc = squeeze(out.fc_FR(:,:,1:length(Gains),1:length(Gg)));
  fc_all(:,:,:,isubj) = reshape(tmp_fc(repmat(mask,[1 1 12 11])),[(76*76-76)/2 12 11]);
  
end



%%

isubj = 26;

figure_w

b=subplot(2,3,1);

outp.corr = outp.corr.*(osc==0); outp.corr(outp.corr==0)=nan;
outp.dist =outp.dist.*(osc==0); outp.dist(outp.dist ==0)=nan;

prc = prctile(reshape(outp.corr_task(isubj,:,:),[1 61*81]),92);
% m=squeeze(outp.corr_task(isubj,:,:)>prc);

par = squeeze(outp.dist(isubj,:,:)); %par(~m)=1.1;
c=imagesc(par,[0.6 1.1]); set(gca,'ydir','normal')
colormap(b,plasma);
axis square tight; set(c,'AlphaData',~isnan(par))

set(gca,'xtick',1:20:81,'xticklabel',{'-0.5';'0';'0.5';'1';'1.5'},'fontsize',8); tp_editplots
set(gca,'ytick',1:20:61,'yticklabel',{'-0.5';'0';'0.5';'1'},'fontsize',8);
line([21 21],[1 21],'color','w','linestyle',':'); line([1 21],[21 21],'color','k','linestyle',':')
line([61 61],[1 41],'color','w','linestyle',':'); line([1 61],[41 41],'color','k','linestyle',':')

xlabel('Change in I'); ylabel('Change in E')


b=subplot(2,3,2);

outp.corr = outp.corr.*(osc==0); outp.corr(outp.corr==0)=nan;
outp.dist =outp.dist.*(osc==0); outp.dist(outp.dist ==0)=nan;

prc = prctile(reshape(outp.corr_task(isubj,:,:),[1 61*81]),92);
m=squeeze(outp.corr_task(isubj,:,:)>prc);

par = squeeze(outp.dist(isubj,:,:)); par(~m)=nan;
c=imagesc(par,[0.6 1.1]); set(gca,'ydir','normal')
colormap(b,plasma);
axis square tight; set(c,'AlphaData',~isnan(par))

set(gca,'xtick',1:20:81,'xticklabel',{'-0.5';'0';'0.5';'1';'1.5'},'fontsize',8); tp_editplots
set(gca,'ytick',1:20:61,'yticklabel',{'-0.5';'0';'0.5';'1'},'fontsize',8);
line([21 21],[1 21],'color','k','linestyle',':'); line([1 21],[21 21],'color','k','linestyle',':')
line([61 61],[1 41],'color','k','linestyle',':'); line([1 61],[41 41],'color','k','linestyle',':')

xlabel('Change in I'); ylabel('Change in E')

a=subplot(2,3,3);

imagesc(fc_task-fc_rest,[-0.03 0.03]);
colormap(a,cmap); axis square off

b=subplot(2,3,4);
imagesc(out.fc_FR(:,:,21,21),[0 0.1])
axis square off; colormap(b,plasma)

b=subplot(2,3,5);
imagesc(out.fc_FR(:,:,41,61),[0 0.1])
axis square off; colormap(b,plasma)

b=subplot(2,3,6);
imagesc(out.fc_FR(:,:,41,61)-out.fc_FR(:,:,21,21),[-0.2 0.2])
axis square off; colormap(b,cmap)


print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_fitting_task_v%d.pdf',v))
