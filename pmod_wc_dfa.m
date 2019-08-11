%% pmod_wc_dfa
% Stochastic simulation of 2*N WC nodes during "rest"
%--------------------------------------------------------------------------

clear

outdir = '~/pmod/proc/dfa/';

% 29-05-2018: fit E and I through resting state recordings. Then obtain
% changes in E and I due to task from recordings and keep those parameters
% fixed for the drug simulations. Vary excitability and gain for the drug
% recordings.
%--------------------------------------------------------------------------
% VERSION 1: 20-10-2018: DETERMINE GLOBAL COUPLING PARAMETER
% %------------------------------------------------------------------------
v           = 1;
Ies         = -14:0.25:2;
Iis         = -14:0.25:2;
Gg          = 1;
Gains       = 0;
nTrials     = 1;
tmax        = 14999;  % in units of tauE
EC          = 0;
%--------------------------------------------------------------------------

% EXCLUDE CERTAIN REGIONS - BCN ordering
k = 1 : 90;
exclude_bcn = [11 15 21 36 37 38 39 52 53 54 55 70 76 80];
include_bcn = find(~ismember(k,exclude_bcn));
% ------------------

% load connectome
if EC
  load ~/pmod/matlab/EC.mat %Matt_EC
  C = EC;
else
  load ~/sc90.mat %Bea SC
  C = SC;
end

C = C/max(C(C>0));
C = C(include_bcn,include_bcn);
N = 76;

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
clear tspan

ds = 10;
Tds = length(0:ds*dt:tmax)-1;
tauEsec = 0.009; % in seconds
resol = ds*dt*tauEsec;
time = (0:ds*dt:tmax-ds*dt)*tauEsec;

sigma = 0.0005;
%Qn = (sigma*dt)^2*eye(2*N);

isub = find( triu(ones(N)) - eye(N) );
%%
for iies = 1: length(Ies)
  iies
  for iiis = 1: length(Iis)
    for iG = 1 : length(Gg)
      %
      fn = sprintf('pmod_wc_dfa_Ie%d_Ii%d_G%d_v%d',iies,iiis,iG,v);
      if tp_parallel(fn,'~/pmod/proc/dfa/',1,0)
        continue
      end
%       
      tic
      g = Gg(iG);
      W = [wEE*eye(N) -wEI*eye(N); wIE*eye(N) -wII*eye(N)];
      
      out.Cee       = zeros(1,1);
      out.CeeSD     = zeros(2,1);
      out.FC_env    = zeros(N,N,1);
      out.Cee_env   = zeros(1,1);
      out.CeeSD_env = zeros(2,1);
      
      out.KOPsd   = zeros(nTrials,1);
      out.KOPmean = zeros(nTrials,1);
      
      % Control params.
      %--------------------
      out.Ie = Ies(iies);
      out.Ii = Iis(iiis);
      %       out.Gain = Gains(igain);
      
      % Working point:
      Io=zeros(2*N,1);
      Io(1:N) = out.Ie;
      Io(N+1:2*N) = out.Ii;
      
      % transfer function:
      gE  = 1;
      gI  = 1;
      aE  = 1/gE;
      aI  = 1/gI;
      Fe  = @(x) 1./(1 + exp(-x/aE) );
      Fi  = @(x) 1./(1 + exp(-x/aI) );
      F   = @(x) [feval(Fe,x(1:N));feval(Fi,x(N+1:2*N))];
      
      T           = Tds*resol; %% define time of interval
      freqs       = (0:Tds/2)/T; %% find the corresponding frequency in Hz
      freq100     = freqs(freqs<100 & freqs>1);
      pp          = 1:10:length(freq100);
      PSD         = zeros(length(pp),N,nTrials);
      frequencies = freq100(1:10:end)';
      
      % RUN SIMULATION
      % ---------------------
      
      for tr = 1:nTrials
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
        z  = rE + 1i*rI;
        
        clear R Ri rI
        
        % KURAMOTO PARAMETERS
        % ---------------------
        %           ku                = sum(z,2)/N;
        %           KOP               = abs(ku);
        %           out.KOPsd(tr,1)   = std(KOP);
        %           out.KOPmean(tr,1) = mean(KOP);
        
        clear ku KOP z
        
        % FC matrix
        % ---------------------
        out.rE(:,tr) = rE;
        %         out.rI(:,tr) = rI;
        out.fs = 1/resol;
        %           out.FC  	= single(out.FC) + single(rc/nTrials);
        flp = 9;           % lowpass frequency of filter
        fhi = 13;
        
        delt  = resol;            % sampling interval
        k     = 4;                  % 2nd order butterworth filter
        fnq   = 1/(2*delt);       % Nyquist frequency
        Wn    = [flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
        [bfilt,afilt] = butter(k,Wn);
        
        % COMPUTE DFA, EXTRACT HURST
        % ---------------------------
        out.env(:,:,tr) = abs(hilbert(filtfilt(bfilt,afilt,rE)));
        
        tmp = tp_dfa(out.env,[3 50],out.fs,0.5,15);
        out.dfa_env = tmp.exp;
        tmp = tp_dfa(rE,[3 50],out.fs,0.5,15);
        out.dfa     = tmp.exp;
        
        clear rE rI env
        toc
      end
      
      save(sprintf('~/pmod/proc/dfa/%s.mat',fn),'out')
      
      % make sure file is saved
      while 1
        try
          load(sprintf('~/pmod/proc/dfa/%s.mat',fn))
          break
        catch me
          save(sprintf('~/pmod/proc/dfa/%s.mat',fn),'out')
        end
      end
      
      tp_parallel(fn,'~/pmod/proc/dfa/',0)
      
    end
  end
end

error('!')
%%

for i = 1:31
  i
    for iG = 1 : 1
      fn = sprintf('pmod_wc_dfa_Ie%d_Ii%d_G%d_v%d',i,i,iG,v);
      
      copyfile(['~/pmod/proc/dfa/' fn '.mat'],'~/pmod/proc/dfa/4guido/')
%       tmp = tp_dfa(out.rE,[3 100],out.fs,0.5,15);
%       dfa(iies,iiis)=mean(tmp.exp);
      
    end
  
end


% v = 1;
%
% ord   = pconn_randomization;
%
% for iies = 1: length(Ies)
%   iies
%   for iiis = 1: length(Iis)
%     for iG = 1 : length(Gg)
%       for igain = 1 : length(Gains)
% %
%
% %         delete(sprintf(['~/pmod/proc/' 'pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d_processing.txt'],iies,iiis,iG,igain,v))
%         if ~exist(sprintf(['~/pmod/proc/' 'pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat'],iies,iiis,iG,igain,v))
%         	delete(sprintf(['~/pmod/proc/' 'pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d_processing.txt'],iies,iiis,iG,igain,v))
%             warning('deleted')
%         end
%
%       end
%     end
%   end
% end
