%% pmod_WC_wholebrain_drugs
% Stochastic simulation of 2*N WC nodes during "rest"
%-------------------------------------------------------------------------

clear

% 29-05-2018: fit E and I through resting state recordings. Then obtain 
% changes in E and I due to task from recordings and keep those parameters
% fixed for the drug simulations. Vary excitability and gain for the drug
% recordings.

%-------------------------------------------------------------------------
% VERSION 1: baseline (no drug), low # trials, long run time (~550s)
%-------------------------------------------------------------------------
% v           = 1;
% Ies         = -4:0.1:-0.8;
% Iis         = -5:0.1:-1;
% Gg          = 0.62; 
% Gains       = 0;
% nTrials     = 1;
% tmax        = 60000; % in units of tauE
%-------------------------------------------------------------------------
% VERSION 2: baseline (no drug), high # trials, short run time (~55s)
%-------------------------------------------------------------------------
v           = 2;
Ies         = -4:0.1:-0.8;
Iis         = -5:0.1:-1;
Gg          = 0.62; 
Gains       = 0;
nTrials     = 20;
tmax        = 10000; % in units of tauE
%-------------------------------------------------------------------------

% load connectome
load ~/pmod/matlab/EC.mat %Matt_EC
C = EC;
C = C/max(C(C>0));
N = size(C,1);

addpath ~/Documents/MATLAB/Colormaps/'Colormaps (5)'/Colormaps/
addpath ~/pconn/matlab
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
%--------------------------------------------------------------------------
% PARAMETER DEFINITIONS
%--------------------------------------------------------------------------

% Connectivity:
wII=4;
wIE=16;
wEI=12;
wEE=12;

tauE = 1;
tauI = 2;
tau = zeros(2*N,1);
tau(1:N) = tauE;
tau(N+1:2*N) = tauI;

dt=0.01;
tspan=0:dt:tmax;
L = length(tspan);

ds = 10;
Tds = length(0:ds*dt:tmax)-1;
tauEsec = 0.009; % in seconds
resol = ds*dt*tauEsec;
time = (0:ds*dt:tmax-ds*dt)*tauEsec;

sigma = 0.0005;
%Qn = (sigma*dt)^2*eye(2*N);

% FILTERS

flp = 8;           % lowpass frequency of filter
fhi = 12;

para.ord = 4;
delt = 1/(1/resol);            % sampling interval
k=4;                  % 2nd order butterworth filter
fnq=1/(2*delt);       % Nyquist frequency
Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt,afilt]=butter(k,Wn);

isub = find( triu(ones(N)) - eye(N) );
%%
for iies = 1:length(Ies)
  for iiis = 1: length(Iis)
    for iG = 1 : length(Gg)
      for igain = 1 : length(Gains)
      
      if ~exist(sprintf(['~/pmod/proc/' 'pmod_WC_wholebrain_drugs_Ie%d_Ii%d_G%d_gain%d_v%d_processing.txt'],iies,iiis,iG,igain,v))
        system(['touch ' '~/pmod/proc/' sprintf('pmod_WC_wholebrain_drugs_Ie%d_Ii%d_G%d_gain%d_v%d_processing.txt',iies,iiis,iG,igain,v)]);
      else
        continue
      end
      tic
      g = Gg(iG);
      W = [wEE*eye(N)+g*C -wEI*eye(N); wIE*eye(N) -wII*eye(N)];
      
      FC = zeros(N,N,1);

      Cee = zeros(1,1);
      CeeSD = zeros(2,1);

      KOPsd   = zeros(nTrials,1);
      KOPmean = zeros(nTrials,1);
      % Control params.
      %--------------------
      Ie = Ies(iies);
      Ii = Iis(iiis);
      out.Ie = Ie;
      out.Ii = Ii;
      out.Gain = Gains(igain);
%       Ie = -3.75;
%       Ii = -4.1;
      
      Io=zeros(2*N,1);
      Io(1:N) = Ie;
      Io(N+1:2*N) = Ii;
      
      % transfer function:
      gE = 1 + Gains(igain);
      gI = 1 + Gains(igain);
      aE  = 1/gE;
      aI  = 1/gI;
      Fe  = @(x) 1./(1 + exp(-x/aE) );
      Fi  = @(x) 1./(1 + exp(-x/aI) );
      F   = @(x) [feval(Fe,x(1:N));feval(Fi,x(N+1:2*N))];
      
      % Working point:
      Io(1:N)     = Ie;
      Io(N+1:2*N) = Ii;
      %
      %     FCval{1}    = zeros(length(isub),nTrials);
      %     Epsilon{1}  = zeros(N,nTrials);
      
      T       = Tds*resol; %% define time of interval
      %     freqs   = (0:Tds/2)/T; %% find the corresponding frequency in Hz
      %     nfreqs  = length(freqs);
      %     freq100 = freqs(freqs<100 & freqs>1);
      %     pp      = 1:10:length(freq100);
      %     PSD     = zeros(length(pp),N,nTrials);
      
      %     PSDstruct(1).frequencies = freq100(1:10:end)';
      
      
      for tr=1:nTrials
        fprintf('Rest, Ie%d, Ii%d, trial%d ...\n',iies,iiis,tr)
        r   = 0.001*rand(2*N,1);
        R   = zeros(Tds,N);
        Ri  = zeros(Tds,N);
        tt  = 0;
        % transient:
        for t = 1:50000
          u = W*r + Io;
          K = feval(F,u);
          r = r + dt*(-r + K)./tau + sqrt(dt)*sigma*randn(2*N,1);
        end
        % simulation
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
        
        %correlation:
        rE = R; 
        rI = Ri;
        
        z             = rE + 1i*rI;
        ku            = sum(z,2)/N;
        KOP           = abs(ku);
        KOPsd(tr,1)   = std(KOP);
        KOPmean(tr,1) = mean(KOP);
        
        rc            = corrcoef(rE);
        
        % COMPUTE LONG RANGE TEMPORAL CORRELATIONS
        % On E time course
        % ---------------------
        tmp             = tp_dfa(R,[3 50],1/resol/ds,0.5,15);
        out.dfa(tr,:,1) = tmp.exp;
        tmp               = tp_dfa(R,[3 50],1/resol,0.5,15);
        out.dfa1(tr,:,1) = tmp.exp;
        % COMPUTE LONG RANGE TEMPORAL CORRELATIONS
        % On envelopes
        % ---------------------
        env = abs(hilbert(filtfilt(bfilt,afilt,rE))); 
%         tmp              = tp_dfa(env,[3 50],1/resol/ds,0.5,15);
%         out.dfa_env(tr,:,1)  = tmp.exp;
        tmp              = tp_dfa(env,[3 50],1/resol,0.5,15);
        out.dfa_env1(tr,:,1) = tmp.exp;
        % ---------------------
        
        out.FC(:,:,1)    	= FC(:,:,1) + rc/nTrials;
        fc           	= rc(isub);
        Cee(1)       	= Cee(1) + mean(fc)/nTrials;
        CeeSD(1)    	= CeeSD(1) + var(fc)/nTrials;
        
        rc              = corrcoef(env);

        out.FC_env(:,:,1)  	= FC(:,:,1) + rc/nTrials;
        fc_env         	= rc(isub);
        out.Cee_env(1)      = Cee(1) + mean(fc)/nTrials;
        out.CeeSD_env(1)    = CeeSD(1) + var(fc)/nTrials;
        %               FCval(:,tr)  = fc;
        
      
        %       rEo       = mean(mean(rE));
        %       rEsd      = var(mean(rE));
        %       Rate(1)   = rEo-rEo;
        %       RateSD(1) = RateSD(1) + rEsd/nTrials;
        
        for i=1:N
          out.osc(tr,:,i) = tp_detect_osc(R(:,i));
          %autocorr
          lags = 1:300;
          acorr = acf(R(:,i),300);
          ii = find(acorr<.2,1,'first');
          if ~isempty(ii)
            out.lags(i,tr) = nan;
          else
          	out.lags(i,tr) = lags(ii);
          end
          %PSD:
%           f = R(:,i) - mean(R(:,i));
%           xdft = fft(f);
%           xdft = xdft(1:floor(Tds/2)+1);
%           pw = (1/(Tds/2)) * abs(xdft).^2;
%           psd = pw(freqs<100 & freqs>1);
%           f = freqs(freqs<100 & freqs>1);
%           fnew = f(1:10:end);
%           psd  = psd(1:10:end);
%           %       pw = gaussfilt(fnew,psd',.5);
%           PSD(:,i,tr) = psd';
        end
        toc
      end

      
      save(sprintf('~/pmod/proc/pmod_WC_wholebrain_drugs_Ie%d_Ii%d_G%d_gain%d_v%d.mat',iies,iiis,iG,igain,v),'out')
      
    end
  end
end
end
error('!')

%% FITTING
% ---------
% LOAD EMPIRICALDFA (LCMV)
% ---------
load(sprintf(['~/pupmod/proc/conn/' 'pupmod_src_dfa_v%d.mat'],1));
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

dfa_emp_rest = nanmean(dfa_all(:,:,1,1,1),2);
% dfa_emp_rest = tp_match_aal(pars,repmat(dfa_emp_rest,[1 90]));
% dfa_emp_rest = dfa_emp_rest(:,1);

dfa_emp_task = nanmean(dfa_all(:,:,1,1,2),2);
% dfa_emp_task = tp_match_aal(pars,repmat(dfa_emp_task,[1 90]));
% dfa_emp_task = dfa_emp_task(:,1);
% ---------

clear dfa r dfa_r dist_fc fc_sim
v=1;
vv =7;
load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));
mask    = logical(tril(ones(90,90),-1));
mask = find(triu(ones(90))-eye(90));

fc_rest =  squeeze(nanmean(cleandat(:,:,:,1,1,6),3));
fc_task =  squeeze(nanmean(cleandat(:,:,:,1,2,6),3));

for iies = 1 : 29%length(Ies)
  iies
  for iiis = 1 : length(Iis)
    for iG = 1
      for igain = 1 : length(Gains)

      %       load(sprintf('~/pmod/proc/pmod_WC_wholebrain_rest_Ie%d_Ii%d_v%d.mat',iies,iiis,vv))
      load(sprintf('~/pmod/proc/pmod_WC_wholebrain_drugs_Ie%d_Ii%d_G%d_gain%d_v%d.mat',iies,iiis,iG,igain,v))
            
      if round(Ies(iies)*100)/100 == -2.8 && round( Iis(iiis)*100)/100 == -3.4
        disp('save stuff')
        FFCC = out.FC;
      end
      
      dfa_sim(:,iies,iiis,iG,igain) = out.dfa_env1(:);
      fc_sim_tmp = tp_match_aal(pars,out.FC(:,:,1));
      
%       r_rest(iies,iiis,iG)=corr(fc_sim_tmp(mask),fc_rest(mask));
      r_rest(iies,iiis,iG,igain) = dot(fc_sim_tmp(mask),fc_rest(mask)) / sqrt(dot(fc_sim_tmp(mask),fc_sim_tmp(mask)) * dot(fc_rest(mask),fc_rest(mask)));
%       r_task(iies,iiis,iG)=corr(fc_sim_tmp(mask),fc_task(mask));
      r_task(iies,iiis,iG,igain) = dot(fc_sim_tmp(mask),fc_task(mask)) / sqrt(dot(fc_sim_tmp(mask),fc_sim_tmp(mask)) * dot(fc_task(mask),fc_task(mask)));

%       kura_dist(iies,iiis,iG)=mean(out.KOPsd)/mean(out.KOPmean)-kura_emp_rest;
      
      idx = find(~isnan(dfa_emp_rest(1:90)'));
      
%       dfa_r_rest (iies,iiis,iG) = dot(dfa_emp_rest(:),dfa_sim(:,iies,iiis,iG)) / sqrt(dot(dfa_sim(:,iies,iiis,iG),dfa_sim(:,iies,iiis,iG)) * dot(dfa_emp_rest(:),dfa_emp_rest(:)))
      dfa_r_rest (iies,iiis,iG,igain) = corr(dfa_emp_rest(:),dfa_sim(:,iies,iiis,iG));
%       dfa_r_rest (iies,iiis,iG) = corr(dfa_emp_rest(:),dfa_sim(:,iies,iiis,iG))*(std(dfa_emp_rest(:))/std(dfa_sim(:,iies,iiis,iG)));
%       dfa_r_rest (iies,iiis,iG) = dot(dfa_emp_rest(:),dfa_sim(:,iies,iiis,iG)) / sqrt(dot(dfa_sim(:,iies,iiis,iG),dfa_sim(:,iies,iiis,iG)) * dot(dfa_emp_task(:),dfa_emp_task(:)))
      dfa_r_task (iies,iiis,iG,igain) = corr(dfa_emp_task(:),dfa_sim(:,iies,iiis,iG));
%       dfa_r_task (iies,iiis,iG) = corr(dfa_emp_task(:),dfa_sim(:,iies,iiis,iG))*(std(dfa_emp_task(:))/std(dfa_sim(:,iies,iiis,iG)));

      dist_fc_rest (iies,iiis,iG,igain)  = mean(fc_sim_tmp(mask))-mean(fc_rest(mask));
      dist_fc_task (iies,iiis,iG,igain)  = mean(fc_sim_tmp(mask))-mean(fc_task(mask));
      
      fc_sim(iies,iiis,iG,igain) = mean(fc_sim_tmp(mask));
      
%       Ies(iies) = out.Ie;
%       Iis(iiis) = out.Ii;
      
      fclose all;

    end
    end
  end
end
% 
% pars.N = 90;
% dfa_emp  =  tp_match_aal(pars,dfa_emp');



%%
idx = [find(round(Ies*100)/100==-2.8) find(round(Iis*100)/100==-3.5000) ];
idx2 = [find(round(Ies*100)/100==-1.8) find(round(Iis*100)/100==-2.4000) ];


% Ie = -2.85;
% Ii = -3.50;
% idx = [-2.25 -4]';
cmap = cbrewer('div', 'RdBu', 100,'pchip');
cmap = cmap(end:-1:1,:);

% i = 4;
h = figure; set(gcf,'color','white'); hold on
set(h,'Position',[10,10,1200,800])
orient(h,'landscape')
gg = 1;

% -------------------------------
% SIMULATED FC MATRIX 
% --------------------------------<<>>

ax1 = subplot(2,4,1);
imagesc(flipud((squeeze(fc_sim(:,:,gg)))),[-0.5  0.5]); 
axis square tight; hold on
scatter(idx(2),length(Ies)-idx(1)+1,20,'markerfacecolor','w','markeredgecolor','k')
scatter(idx2(2),length(Ies)-idx2(1)+1,20,'markerfacecolor','r','markeredgecolor','k')

title('FC_{sim}')
colormap(ax1,cmap)
tp_editplots;

% -------------------------------
% SIMULATED CORRELATIONS
% -------------------------------
ax1 = subplot(2,4,2);
imagesc(flipud((squeeze(fc_sim(:,:,gg))-squeeze(fc_sim(:,:,gg)))),[-0.01  0.01]); 
axis square tight; hold on
scatter(idx(2),length(Ies)-idx(1)+1,20,'markerfacecolor','w','markeredgecolor','k')
scatter(idx2(2),length(Ies)-idx2(1)+1,20,'markerfacecolor','r','markeredgecolor','k')

title('FC_{sim}')
colormap(ax1,cmap)
tp_editplots;


ax5=subplot(2,4,5);

% mask3 = squeeze(mean(dfa_sim(:,:,:,gg)))>0.5 & squeeze(mean(dfa_sim(:,:,:,gg)))<0.9;
imagesc(flipud((squeeze(mean(dfa_sim(:,:,:,gg)))-squeeze(mean(dfa_sim(:,:,:,gg))))),[-0.01 0.01]);  
axis square tight; hold on
scatter(idx(2),length(Ies)-idx(1)+1,20,'markerfacecolor','w','markeredgecolor','k')
scatter(idx2(2),length(Ies)-idx2(1)+1,20,'markerfacecolor','r','markeredgecolor','k')


colormap(ax5,cmap)

title('\alpha_{sim}')
tp_editplots
hold on

% -------------------------------
% CORRELATION SIMULATED DFA / RESTING DFA
% -------------------------------

ax4=subplot(2,4,6);
mask4 = squeeze(dfa_r_rest(:,:,gg))>-1;
imagesc(flipud(squeeze(dfa_r_rest(:,:,gg))),[-0.3 0.3]);  
axis square tight; hold on

colormap(ax4,cmap)
title('Correlation (\alpha_{sim},\alpha_{emp})')
tp_editplots; 

scatter(idx(2),length(Ies)-idx(1)+1,20,'markerfacecolor','w','markeredgecolor','k')
scatter(idx2(2),length(Ies)-idx2(1)+1,20,'markerfacecolor','r','markeredgecolor','k')

% -------------------------------
% DISTANCE SIMULATED DFA - RESTING DFA
% -------------------------------

ax5=subplot(2,4,7);

mask5 = abs(squeeze(mean(dfa_sim(:,:,:,gg),1))-mean(dfa_emp_rest))<0.5;
imagesc(flipud((squeeze(mean(dfa_sim(:,:,:,gg),1))-mean(dfa_emp_rest))),[-0.05 0.05]);  
axis square tight; hold on
colormap(ax5,cmap)

title('Difference (\alpha_{Sim} - \alpha_{Emp})')
tp_editplots; hold on
scatter(idx(2),length(Ies)-idx(1)+1,20,'markerfacecolor','w','markeredgecolor','k')
scatter(idx2(2),length(Ies)-idx2(1)+1,20,'markerfacecolor','r','markeredgecolor','k')


% -------------------------------
% MASKS RESTING DFA
% -------------------------------
% 
% ax5=subplot(2,4,8);
% imagesc(flipud(mask3&mask4&mask5),[-1 1]);  
% axis square tight;hold on
% scatter(idx(2),length(Ies)-idx(1)+1,20,'markerfacecolor','w','markeredgecolor','k')
% % scatter(idx2(2),length(Ies)-idx2(1)+1,20,'markerfacecolor','r','markeredgecolor','k')
% 
% colormap(ax5,cmap)
% tp_editplots

% -------------------------------
% ADD AXIS LABELS TO ALL SUBPLOTS
% -------------------------------

ax = findobj(h,'type','axes');
for iax = 1 : length(ax)
  if length(Iis)== 1
    set(ax(iax),'xtick',1:4:length(Ies),'xticklabel',Ies(1:4:length(Ies)),'tickdir','out')
    set(ax(iax),'ytick',1:3:length(Gg),'yticklabel',Gg(length(Gg):-3:1),'tickdir','out')
    xlabel(ax(iax),'Exc. input I_E'); ylabel(ax(iax),'Global coupling G');
  else
    set(ax(iax),'xtick',1:round(length(Iis)/5):length(Iis),'xticklabel',Iis(1:round(length(Iis)/5):length(Iis)),'tickdir','out')
    set(ax(iax),'ytick',1:3:length(Ies),'yticklabel',Ies(length(Ies):-3:1),'tickdir','out')
    xlabel(ax(iax),'Inh. input I_I'); ylabel(ax(iax),'Exc. input I_E'); 
  end
end

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_WC_wholebrain_rest_drugs_sub_v%d.pdf',v))

title('Masked values')
tp_editplots
hold on

% BAR PLOTS
%%
figure; set(gcf,'color','w');

diff_atx(1) = fc_sim(idx(1),idx(2),1,3)-fc_sim(idx(1),idx(2),1,2);
diff_atx(1) = fc_sim(idx2(1),idx2(2),1,1)-fc_sim(idx2(1),idx2(2),1,2);

diff_dpz(1) = fc_sim(idx(1),idx(2),1,3)-fc_sim(idx(1),idx(2),1,2);
diff_dpz(2) = fc_sim(idx2(1),idx2(2),1,1)-fc_sim(idx2(1),idx2(2),1,2);

bar([1 2 3 4],[diff_atx diff_dpz])

%%
% -------------------------------
% COMBINED MASKS
% -------------------------------

figure;
imagesc(flipud(mask1&mask2&mask3&mask4&mask5),[-1 1]);  axis square tight
colormap(cmap)

if length(Iis)== 1
  set(gca,'xtick',1:4:length(Ies),'xticklabel',Ies(1:4:length(Ies)),'tickdir','out')
  set(gca,'ytick',1:3:length(Gg),'yticklabel',Gg(length(Gg):-3:1),'tickdir','out')
  xlabel('Exc. input I_E'); ylabel('Global coupling G');
else
  set(gca,'xtick',1:round(length(Iis)/5):length(Iis),'xticklabel',Iis(1:round(length(Iis)/5):length(Iis)),'tickdir','out')
  set(gca,'ytick',1:3:length(Ies),'yticklabel',Ies(length(Ies):-3:1),'tickdir','out')
  xlabel('Inh. input I_I'); ylabel('Exc. input I_E'); 
end
hold on
tp_editplots
scatter(idx(2),length(Ies)-idx(1)+1,20,'markerfacecolor','w','markeredgecolor','k')
scatter(idx2(2),length(Ies)-idx2(1)+1,20,'markerfacecolor','r','markeredgecolor','k')


print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_WC_wholebrain_rest_mask_v%d.pdf',v))

% --------------------------------------------------------------------


%% GET TASK PARAMETERS
%
% idx = [find(Ies==-2.9) find(Iis==-5.2) ];
% idx = [-2.25 -4]';
cmap = cbrewer('div', 'RdBu', 100,'pchip');
cmap = cmap(end:-1:1,:);

% i = 4;
h = figure; set(gcf,'color','white'); hold on
set(h,'Position',[10,10,1200,800])
orient(h,'landscape')

% -------------------------------
% SIMULATED FC MATRIX 
% -------------------------------

ax1 = subplot(2,4,1);
mask1 = fc_sim(:,:,gg) > 0.10 & fc_sim(:,:,gg) < 0.8;
imagesc(flipud((squeeze(fc_sim(:,:,gg)))),[0  1]); 
axis square tight; hold on
scatter(idx(2),length(Ies)-idx(1)+1,20,'markerfacecolor','w','markeredgecolor','k')

title('FC_{sim}')
colormap(ax1,inferno)
tp_editplots;

% -------------------------------
% CORRELATION OF SIMULATED FC WITH RESTING FC 
% -------------------------------

ax2 = subplot(2,4,2);
mask2 = r_task(:,:,gg)>0.2;
imagesc(flipud(squeeze(r_task(:,:,gg))),[-0.3 0.3]); 
axis square tight; hold on
scatter(idx(2),length(Ies)-idx(1)+1,20,'markerfacecolor','w','markeredgecolor','k')
colormap(ax2,cmap)
title('Correlation (FC_{sim},FC_{emp})'); hold on
tp_editplots; 

% -------------------------------
% DISTANCE OF SIMULATED FC FROM RESTING FC 
% -------------------------------

ax3=subplot(2,4,3); 

imagesc(flipud(squeeze(dist_fc_task(:,:,gg))),[-0.5 0.5]);  
axis square tight; hold on

scatter(idx(2),length(Ies)-idx(1)+1,20,'markerfacecolor','w','markeredgecolor','k')
colormap(ax3,cmap)
title('Difference (FC_{sim} - FC_{emp})')
tp_editplots;

% -------------------------------
% MASK: FC
% -------------------------------

ax4=subplot(2,4,4);

imagesc(flipud((mask1&mask2)),[-1 1]);  
axis square tight; hold on
scatter(idx(2),length(Ies)-idx(1)+1,20,'markerfacecolor','w','markeredgecolor','k')
colormap(ax4,cmap)
title('Masked values')
tp_editplots; 

% -------------------------------
% SIMULATED SCALING EXPONENT
% -------------------------------

ax5=subplot(2,4,5);

mask3 = squeeze(mean(dfa_sim(:,:,:,gg)))>0.5 & squeeze(mean(dfa_sim(:,:,:,gg)))<1;
imagesc(flipud((squeeze(mean(dfa_sim(:,:,:,gg))))),[0.5 1]);  
axis square tight; hold on
scatter(idx(2),length(Ies)-idx(1)+1,20,'markerfacecolor','w','markeredgecolor','k')

colormap(ax5,inferno)

title('\alpha_{sim}')
tp_editplots
hold on

% -------------------------------
% CORRELATION SIMULATED DFA / RESTING DFA
% -------------------------------

ax4=subplot(2,4,6);
mask4 = squeeze(dfa_r_task(:,:,gg))>0.1;
imagesc(flipud(squeeze(dfa_r_task(:,:,gg))),[-0.5 0.5]);  
axis square tight; hold on

colormap(ax4,cmap)
title('Correlation (\alpha_{sim},\alpha_{emp})')
tp_editplots; 

scatter(idx(2),length(Ies)-idx(1)+1,20,'markerfacecolor','w','markeredgecolor','k')

% -------------------------------
% DISTANCE SIMULATED DFA - RESTING DFA
% -------------------------------

ax5=subplot(2,4,7);

mask5 = abs(squeeze(mean(dfa_sim(:,:,:,gg),1))-mean(dfa_emp_task))<0.1;
imagesc(flipud((squeeze(mean(dfa_sim(:,:,:,gg),1))-mean(dfa_emp_task))),[-0.2 0.2]);  
axis square tight; hold on
colormap(ax5,cmap)

title('Difference (\alpha_{Sim} - \alpha_{Emp})')
tp_editplots; hold on
scatter(idx(2),length(Ies)-idx(1)+1,20,'markerfacecolor','w','markeredgecolor','k')

% -------------------------------
% MASKS RESTING DFA
% -------------------------------

ax5=subplot(2,4,8);
imagesc(flipud(mask3&mask4&mask5),[-1 1]);  
axis square tight;hold on
scatter(idx(2),length(Ies)-idx(1)+1,20,'markerfacecolor','w','markeredgecolor','k')
colormap(ax5,cmap)
tp_editplots

% -------------------------------
% ADD AXIS LABELS TO ALL SUBPLOTS
% -------------------------------

ax = findobj(h,'type','axes');
for iax = 1 : length(ax)
  if length(Iis)== 1
    set(ax(iax),'xtick',1:4:length(Ies),'xticklabel',Ies(1:4:length(Ies)),'tickdir','out')
    set(ax(iax),'ytick',1:3:length(Gg),'yticklabel',Gg(length(Gg):-3:1),'tickdir','out')
    xlabel(ax(iax),'Exc. input I_E'); ylabel(ax(iax),'Global coupling G');
  else
    set(ax(iax),'xtick',1:round(length(Iis)/5):length(Iis),'xticklabel',Iis(1:round(length(Iis)/5):length(Iis)),'tickdir','out')
    set(ax(iax),'ytick',1:3:length(Ies),'yticklabel',Ies(length(Ies):-3:1),'tickdir','out')
    xlabel(ax(iax),'Inh. input I_I'); ylabel(ax(iax),'Exc. input I_E'); 
  end
end

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_WC_wholebrain_task_sub_v%d.pdf',v))

title('Masked values')
tp_editplots
hold on
%%
% -------------------------------
% COMBINED MASKS
% -------------------------------

figure;
imagesc(flipud(mask1&mask2&mask3&mask4&mask5),[-1 1]);  axis square tight
colormap(cmap)

if length(Iis)== 1
  set(gca,'xtick',1:4:length(Ies),'xticklabel',Ies(1:4:length(Ies)),'tickdir','out')
  set(gca,'ytick',1:3:length(Gg),'yticklabel',Gg(length(Gg):-3:1),'tickdir','out')
  xlabel('Exc. input I_E'); ylabel('Global coupling G');
else
  set(gca,'xtick',1:round(length(Iis)/5):length(Iis),'xticklabel',Iis(1:round(length(Iis)/5):length(Iis)),'tickdir','out')
  set(gca,'ytick',1:3:length(Ies),'yticklabel',Ies(length(Ies):-3:1),'tickdir','out')
  xlabel('Inh. input I_I'); ylabel('Exc. input I_E'); 
end
hold on
tp_editplots
scatter(idx(2),length(Ies)-idx(1)+1,20,'markerfacecolor','w','markeredgecolor','k')

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_WC_wholebrain_task_mask_v%d.pdf',v))

% --------------------------------------------------------------------


%% TESTS WITH KURAMOTO
% 
% figure;
% 
% imagesc(flipud(kura_dist(:,:,gg)),[-0.5 0.5]);  
% axis square tight; hold on
% scatter(idx(2),length(Ies)-idx(1)+1,20,'markerfacecolor','w','markeredgecolor','k')
% colormap(cmap)
% title('Masked values')
% tp_editplots; 


pars = [];
pars.grid = 'medium';
pars.N = 90;

dfa_emp_rest = tp_match_aal(pars,fc_rest);


