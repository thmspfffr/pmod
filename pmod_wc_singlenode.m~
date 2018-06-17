%% pmod_wc_singlenode
% Stochastic simulation of 2*N WC nodes during "rest"
%-------------------------------------------------------------------------

clear

%-------------------------------------------------------------------------
% VERSION 1
%-------------------------------------------------------------------------
v = 1;
Ies = -8:0.1:2;
Iis = -8:0.1:2;
Gg  = 0; % 0.62
nTrials = 1;
envelopes = 1;
N = 5;%size(C,1);
tmax  = 100000; % in units of tauE
%-------------------------------------------------------------------------
% VERSION 2
%-------------------------------------------------------------------------
v = 2;
Ies = -8:0.1:2;
Iis = -8:0.1:2;
Gg  = 0; % 0.62
nTrials = 1;
envelopes = 1;
N = 1;%size(C,1);
tmax  = 1000000; % in units of tauE
%-------------------------------------------------------------------------

%-----
% load connectome
load ~/pmod/matlab/EC.mat %Matt_EC
C = EC;
C = C/max(C(C>0));

C=C(1:N,1:N);

addpath ~/Documents/MATLAB/Colormaps/'Colormaps (5)'/Colormaps/
addpath ~/pconn/matlab
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
%--------------------------------------------------------------------------
% PARAMETER DEFINITIONS
%--------------------------------------------------------------------------
% Connectivity:
wII = 4;
wIE = 16;
wEI = 12;
wEE = 12;

tauE          = 1;
tauI          = 2;
tau           = zeros(2*N,1);
tau(1:N)      = tauE;
tau(N+1:2*N)  = tauI;

dt    = 0.01;
tspan = 0:dt:tmax;
L     = length(tspan);

ds      = 10;
Tds     = length(0:ds*dt:tmax)-1;
tauEsec = 0.009; % in seconds
resol   = ds*dt*tauEsec;
time    = (0:ds*dt:tmax-ds*dt)*tauEsec;

sigma = 0.0005;

% FILTERS

flp = 8;           % lowpass frequency of filter
fhi = 12;

para.ord = 4;
delt = 1/(1/resol/ds);            % sampling interval
k=4;                  % 2nd order butterworth filter
fnq=1/(2*delt);       % Nyquist frequency
Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt,afilt]=butter(k,Wn);


% dI_drug(N+1:2*N) = dIi_drug;

isub = find( triu(ones(N)) - eye(N) );
%%
for iies = 1:length(Ies)
  for iiis = 1: length(Iis)
    for iG = 1 : length(Gg)
% % % %       
      if ~exist(sprintf(['~/pmod/proc/' 'pmod_wc_singlenode_Ie%d_Ii%d_G%d_v%d_processing.txt'],iies,iiis,iG,v))
        system(['touch ' '~/pmod/proc/' sprintf('pmod_wc_singlenode_Ie%d_Ii%d_G%d_v%d_processing.txt',iies,iiis,iG,v)]);
      else
        continue
      end
%       Ie = 0;
%       Ii = -5;
      g = Gg(iG);
      W = [wEE*eye(N)+g*C -wEI*eye(N); wIE*eye(N) -wII*eye(N)];
      
      % Control params.
      %--------------------
      Ie = Ies(iies);
      Ii = Iis(iiis);
                                                                                                                                             
      Io=zeros(2*N,1);
%       Io(1:N) = Ie;
%       Io(N+1:2*N) = Ii;
      
      % transfer function:
      gE  = 1;
      gI  = 1;
      aE  = 1/gE;
      aI  = 1/gI;
      Fe  = @(x) 1./(1 + exp(-x/aE) );
      Fi  = @(x) 1./(1 + exp(-x/aI) );
      F   = @(x) [feval(Fe,x(1:N));feval(Fi,x(N+1:2*N))];
      
      % Working point:
      Io(1:N)     = Ie;
      Io(N+1:2*N) = Ii;

      T       = Tds*resol; %% define time of interval

      for tr=1:nTrials
        fprintf('Rest, Ie%d, Ii%d, trial%d ...\n',iies,iiis,tr)
        r   = 0.001*rand(2*N,1);
        R   = zeros(Tds,N);
        Ri  = zeros(Tds,N);
        tt  = 0;
        % transient:
        for t = 1:25000
          u = W*r + Io;
          K = feval(F,u);
          r = r + dt*(-r + K)./tau + sqrt(dt)*sigma*randn(2*N,1);
        end
        % simulation
        for t = 1:L
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
        rE = R; clear R
        rI = Ri;clear Ri
%         plot(rE(:,1)); hold on
%         drawnow
      	out.osc(tr) = tp_detect_osc(rE(:,1));
        env = abs(hilbert(filtfilt(bfilt,afilt,rE))); 
%         env = resample(env,1,50);

        tmp  = tp_dfa(env,[3 50],1/resol,0.5,15);
        out.dfa_env(:,tr) = tmp.exp;

%          rE = resample(rE,1,50);
        tmp 	= tp_dfa(rE,[3 50],1/resol,0.5,15);
        out.dfa_rE(:,tr) = tmp.exp;
   
        
        out.freq = 0:(1/resol)/length(rE):(1/resol)/2;
        p = abs(fft(rE)); p = p (1:size(p,1)/2+1,:);
        out.pow = p; clear p
        
        p = abs(fft(env)); p = p (1:size(p,1)/2+1,:);
        out.envpow = p;
        
        out.freq2 = 0:(1/resol/ds)/length(rE):(1/resol/ds)/2;
        p = abs(fft(rE)); p = p (1:size(p,1)/2+1,:);
        out.pow2 = p; clear p
        
        p = abs(fft(env)); p = p (1:size(p,1)/2+1,:);
        out.envpow2 = p;

      end
      
      save(sprintf('~/pmod/proc/pmod_wc_singlenode_Ie%d_Ii%d_G%d_v%d.mat',iies,iiis,iG,v),'out')
      
    end
  end
end

% error('!')

%%
if ~exist(sprintf('~/pmod/proc/pmod_wc_singlenode_alloutp_v%d.mat',v))

  for iies = 1: length(Ies)
    iies
    for iiis = 1: length(Iis)
      iiis
      for iG = 1 : length(Gg)
     try
       load(sprintf('~/pmod/proc/pmod_wc_singlenode_Ie%d_Ii%d_G%d_v%d.mat',iies,iiis,iG,v))
     catch me
       delete(sprintf('~/pmod/proc/pmod_wc_singlenode_Ie%d_Ii%d_G%d_v%d_processing.txt',iies,iiis,iG,v))
     end
       allout.dfa1(iies,iiis) = mean(out.dfa_rE);
       allout.dfa2(iies,iiis) = mean(out.dfa_rE2);
       allout.dfa_env1(iies,iiis) = mean(out.dfa_env);
       allout.dfa_env2(iies,iiis) = mean(out.dfa_env2);
       allout.osc(iies,iiis) = out.osc;
  %      dfa_z(iies,iiis) = mean(out.dfa_z);
      end
    end
  end

  save(sprintf('~/pmod/proc/pmod_wc_singlenode_alloutp_v%d.mat',v),'allout')

else 
  
	load(sprintf('~/pmod/proc/pmod_wc_singlenode_alloutp_v%d.mat',v))
end

%%
addpath ~/Documents/MATLAB/Colormaps/Colormaps' (5)'/Colormaps/
clim = [0.4 0.7];
dfa=allout.dfa1;

mask = dfa>clim(1) & dfa<clim(2);


h=figure; subplot(1,2,1);
ft_plot_matrix(allout.osc,'clim',clim,'highlight',mask,'highlightstyle','opacity')
colormap(plasma)
set(gca,'YDir','normal')
set(gca,'xtick',1:10:length(Iis),'xticklabel',[Iis(1:10:length(Iis))],'tickdir','out')
set(gca,'ytick',1:10:length(Ies),'yticklabel',[Iis(1:10:length(Ies))],'tickdir','out')
xlabel('Inhibitory input I_i'); ylabel('Excitatory input I_E'); 
axis square
colorbar
subplot(1,2,2);
dfa=allout.dfa2;

clim = [0.4 0.7];
mask = dfa>clim(1) & dfa<clim(2);

ft_plot_matrix(dfa,'clim',clim,'highlight',mask,'highlightstyle','opacity')
colormap(plasma)
set(gca,'YDir','normal')
set(gca,'xtick',1:10:length(Iis),'xticklabel',[Iis(1:10:length(Iis))],'tickdir','out')
set(gca,'ytick',1:10:length(Ies),'yticklabel',[Iis(1:10:length(Ies))],'tickdir','out')
% set(gca,'YDir','normal')
xlabel('Inhibitory input I_i'); ylabel('Excitatory input I_E'); 

axis square
colorbar
print(gcf,'-dpdf',sprintf('~/fc_sim.pdf'))
% load(sprintf(['~/pupmod/proc/conn/' 'pupmod_src_dfa_v%d.mat'],1));
% dfa_rest = squeeze(nanmean(dfa_all(:,:,1,1,1),2));
% dfa_task = squeeze(nanmean(dfa_all(:,:,1,1,2),2));

%%
% figure; set(gcf,'color','w')
% imagesc((dfa-dfa_rest(1)),[-0.11 0.11])
% colormap(cmap)
% set(gca,'YDir','normal')
% set(gca,'xtick',1:10:length(Iis),'xticklabel',[Iis(1:10:length(Iis))],'tickdir','out')
% set(gca,'ytick',1:10:length(Ies),'yticklabel',[Iis(1:10:length(Ies))],'tickdir','out')
% xlabel('Inhibition'); ylabel('Excitation')
figure; set(gcf,'color','w'); hold on
imagesc((dfa-dfa_rest(1)),[-0.10 0.10])
colormap(cmap)
set(gca,'YDir','normal')
set(gca,'xtick',1:10:length(Iis),'xticklabel',[Iis(1:10:length(Iis))],'tickdir','out')
set(gca,'ytick',1:10:length(Ies),'yticklabel',[Ies(1:10:length(Ies))],'tickdir','out')
xlabel('Inhibition'); ylabel('Excitation')

% % set(gca,'YDir','normal')
for i = 1 : 90
[minMatrix]=min(abs(dfa(:)-dfa_rest(i)))
[row, col] = find(dfa-dfa_rest(i)==minMatrix);
scatter(col,row,20,'markerfacecolor','w','markeredgecolor','k')
end

for i = 1 : 90
[minMatrix]=min(abs(dfa(:)-dfa_task(i)))
[row, col] = find(dfa-dfa_task(i)==minMatrix);
scatter(col,row,20,'markerfacecolor','r','markeredgecolor','k')
end

%% 
