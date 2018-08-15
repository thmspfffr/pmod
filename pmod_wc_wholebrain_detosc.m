%% pmod_wc_wholebrain_detosc
% Determine oscillatory regime
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
% Ies         = -4:0.1:-1;
% Iis         = -5:0.1:-1;
% Gg          = 0.62;
% Gains       = -0.5:0.25:0.5;
% nTrials     = 3;
% tmax        = 5000; % in units of tauE
% N           = 90;
% wII=4;
% wIE=16;
% wEI=12;
% wEE=12;
%-------------------------------------------------------------------------
% VERSION 2: woolrich 1
%-------------------------------------------------------------------------
v           = 2;
Ies         = -10:0.5:10;
Iis         = -10:0.5:10;
Gg          = 0:0.1:1;
Gains       = 0;
nTrials     = 1;
tmax        = 1000; % in units of tauE
wins        = [3 50];
N           = 90;
% wII=4;
% wIE=16;
% wEI=12;
% wEE=3.5;
%-------------------------------------------------------------------------
% VERSION 3:  woolrich 2
%-------------------------------------------------------------------------
v           = 3;
Ies         = -7:0.5:8;
Iis         = -10:0.5:4;
Gg          = 0.3:0.2:0.9;
Gains       = 0;
nTrials     = 1;
tmax        = 10000; % in units of tauE
wins        = [3 50];
wII=4;
wIE=16;
wEI=12;
wEE=12;
%-------------------------------------------------------------------------

% load connectome
load ~/pmod/matlab/EC.mat %Matt_EC
C = EC;
C = C/max(C(C>0));
if N == 1
  C = 0;
end
% N = size(C,1);

addpath ~/Documents/MATLAB/Colormaps/'Colormaps (5)'/Colormaps/
addpath ~/pconn/matlab
addpath ~/Documents/MATLAB/cbrewer/cbrewer/
%--------------------------------------------------------------------------
% PARAMETER DEFINITIONS
%--------------------------------------------------------------------------

% Connectivity:


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
for iies = 1: length(Ies)
  for iiis = 1:length(Iis)
    for iG = 1 : length(Gg)
      for igain = 1 : length(Gains)
%         
        if ~exist(sprintf(['~/pmod/proc/' 'pmod_wc_wholebrain_detosc_Ie%d_Ii%d_G%d_gain%d_v%d_processing.txt'],iies,iiis,iG,igain,v))
          system(['touch ' '~/pmod/proc/' sprintf('pmod_wc_wholebrain_detosc_Ie%d_Ii%d_G%d_gain%d_v%d_processing.txt',iies,iiis,iG,igain,v)]);
        else
          continue
        end
        tic
        g = Gg(iG);
        W = [wEE*eye(N)+g*C -wEI*eye(N); wIE*eye(N) -wII*eye(N)];
        
%         FC = zeros(N,N,1);
        
%         Cee = zeros(1,1);
%         CeeSD = zeros(2,1);
        
%         KOPsd   = zeros(nTrials,1);
%         KOPmean = zeros(nTrials,1);
        
        %--------------------
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
        
%         T       = Tds*resol; %% define time of interval
%         freqs   = (0:Tds/2)/T; %% find the corresponding frequency in Hz
% %         nfreqs  = length(freqs);
%         freq100 = freqs(freqs<100 & freqs>1);
%         pp      = 1:10:length(freq100);
%         PSD     = zeros(length(pp),N,nTrials);
        
%         PSDstruct(1).frequencies = freq100(1:10:end)';
        
        
        for tr=1:nTrials
          fprintf('Rest, Ie%d, Ii%d, trial%d ...\n',iies,iiis,tr)
          r   = 0.001*rand(2*N,1);
          R   = zeros(Tds,N);
          Ri  = zeros(Tds,N);
          tt  = 0;
          %         transient:
          for t = 1:500
            u = W*r + Io;
            K = feval(F,u);
            r = r + dt*(-r + K)./tau; %+ sqrt(dt);
          end
          %         simulation
          for t = 1:L
            %           100*t/L
            u = W*r + Io;
            K = feval(F,u);
            r = r + dt*(-r+ K)./tau; %+ sqrt(dt);
            if mod(t,ds)==0
              tt=tt+1;
              R(tt,:)  = r(1:N);
              Ri(tt,:)  = r(N+1:end);
            end
          end
          
          %correlation:
%           rE = R;
%           rI = Ri;
          
%           plot(rE(:,1));  hold on
          
%           z             = rE + 1i*rI;
          %         ku            = sum(z,2)/N;
          %         KOP           = abs(ku);
          %         KOPsd(tr,1)   = std(KOP);
          %         KOPmean(tr,1) = mean(KOP);
          
          %         rc            = corrcoef(rE);
          % ---------------------
          % COMPUTE LONG RANGE TEMPORAL CORRELATIONS
          % On E time course
          % ---------------------
          %         tmp               = tp_dfa(R,[5 50],1/resol,0.5,15);
          %         out.dfa1(tr,:,1) = tmp.exp;
          % ---------------------
          % COMPUTE LONG RANGE TEMPORAL CORRELATIONS
          % On envelopes
          % ---------------------
          %         env = abs(hilbert(filtfilt(bfilt,afilt,rE)));
          %         tmp              = tp_dfa(env,[3 50],1/resol/ds,0.5,15);
          %         out.dfa_env(tr,:,1)  = tmp.exp;
          %         tmp              = tp_dfa(env,[5 50],1/resol,0.5,15);
          %         out.dfa_env1(tr,:,1) = tmp.exp;
          % ---------------------
          %
          %         out.FC(:,:,1)    	= FC(:,:,1) + rc/nTrials;
          %         fc           	= rc(isub);
          %         Cee(1)       	= Cee(1) + mean(fc)/nTrials;
          %         CeeSD(1)    	= CeeSD(1) + var(fc)/nTrials;
          %
          %         rc              = corrcoef(env);
          %
          %         out.FC_env(:,:,1)  	= FC(:,:,1) + rc/nTrials;
          %         fc_env         	= rc(isub);
          %         out.Cee_env(1)      = Cee(1) + mean(fc)/nTrials;
          %         out.CeeSD_env(1)    = CeeSD(1) + var(fc)/nTrials;
          %               FCval(:,tr)  = fc;
          
          %       rEo       = mean(mean(rE));
          %       rEsd      = var(mean(rE));
          %       Rate(1)   = rEo-rEo;
          %       RateSD(1) = RateSD(1) + rEsd/nTrials;
          
          for i=1:N
            [out.osc1(tr,:,i)  ] = tp_detect_osc(R(:,i));    
          end
          toc
        end
        
        
        save(sprintf('~/pmod/proc/pmod_wc_wholebrain_detosc_Ie%d_Ii%d_G%d_gain%d_v%d.mat',iies,iiis,iG,igain,v),'out')
        
      end
    end
  end
end
error('!')

%%
osc1 = zeros(length(Ies),length(Iis),length(Gg),length(Gains));
vv =3;
for iies = 1 : length(Ies)
  iies
  for iiis = 1 : length(Iis)
    for iG = 1:length(Gg)
      for igain = 1%:length(Gains)
        %       load(sprintf('~/pmod/proc/pmod_WC_wholebrain_rest_Ie%d_Ii%d_v%d.mat',iies,iiis,vv))
        load(sprintf('~/pmod/proc/pmod_wc_wholebrain_detosc_Ie%d_Ii%d_G%d_gain%d_v%d.mat',iies,iiis,iG,igain,vv))
        
        osc1(iies,iiis,iG,igain) = mean(squeeze(mean(squeeze(out.osc1),1)));
%             osc2(iies,iiis,iG,igain) = mean(squeeze(mean(squeeze(out.osc2),1)));
%         osc3(iies,iiis,iG,igain) = mean(squeeze(mean(squeeze(out.osc3),1)));
%         catch me
%           osc1(iies,iiis,iG,igain) = nan;
%         end
      end
    end
  end
end