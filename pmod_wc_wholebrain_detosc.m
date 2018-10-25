%% pmod_wc_wholebrain_detosc
% Determine oscillatory regime
%-------------------------------------------------------------------------

clear

% 29-05-2018: fit E and I through resting state recordings. Then obtain
% changes in E and I due to task from recordings and keep those parameters
% fixed for the drug simulations. Vary excitability and gain for the drug
% recordings.

%-------------------------------------------------------------------------
% VERSION 1: 20-10-2018
% %-------------------------------------------------------------------------
% v           = 1;
% Ies         = -4:0.2:-1;
% Iis         = -5:0.2:-2;
% Gg          = 0:0.01:1.2;
% Gains       = 0; % 0:0.05:0.2, 3 trials, 
% nTrials     = 1;
% tmax        = 5000; % in units of tauE
%-------------------------------------------------------------------------
% VERSION 2: 20-10-2018
%-------------------------------------------------------------------------
v           = 2;
Ies         = -4:0.025:-1;
Iis         = -5:0.025:-2;
Gg          = 0.87; % this is where correlation peaks 
Gains       = 0; 
nTrials     = 1;
tmax        = 6500; % in units of tauE
%-------------------------------------------------------------------------

% load connectome
load ~/pmod/matlab/EC.mat %Matt_EC
% EXCLUDE CERTAIN REGIONS - BCN ordering
k = 1 : 90;
exclude_bcn = [11 15 21 36 37 38 39 52 53 54 55 70 76 80];
include_bcn = find(~ismember(k,exclude_bcn));
C = EC(include_bcn,include_bcn);
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
        
        fn = sprintf('pmod_wc_wholebrain_detosc_Ie%d_Ii%d_G%d_gain%d_v%d',iies,iiis,iG,igain,v);
        if tp_parallel(fn,'~/pmod/proc/',1)
          continue
        end

        tic
        g = Gg(iG);
        W = [wEE*eye(N)+g*C -wEI*eye(N); wIE*eye(N) -wII*eye(N)];

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

          for i=1:N
            [out.osc1(tr,:,i)  ] = tp_detect_osc(R(:,i));    
          end
          toc
        end
               
        save(sprintf('~/pmod/proc/%s.mat',fn),'out')
        
        while 1
          try
            load(sprintf('~/pmod/proc/%s.mat',fn))
          catch me
            save(sprintf('~/pmod/proc/%s.mat',fn),'out')
            continue
          end
          break
        end
  
        tp_parallel(fn,'~/pmod/proc/',0)

      end
    end
  end
end
error('!')

%% LOAD / DELETE ALL FILES
vv = 2;

if ~exist(sprintf('~/pmod/proc/pmod_wc_wholebrain_detosc_all_v%d.mat',vv))
  if vv ==4
    osc1 = zeros(length(Ies),length(Iis),length(Gg),1);
  else
    osc1 = zeros(length(Ies),length(Iis),length(Gg),length(Gains));
  end
  for iies = 1 : length(Ies)
    iies
    for iiis = 1 : length(Iis)
      for iG = 1:length(Gg)
        for igain = 1:length(Gains)
          
          %       load(sprintf('~/pmod/proc/pmod_WC_wholebrain_rest_Ie%d_Ii%d_v%d.mat',iies,iiis,vv))
          load(sprintf('~/pmod/proc/pmod_wc_wholebrain_detosc_Ie%d_Ii%d_G%d_gain%d_v%d.mat',iies,iiis,iG,igain,vv))
          
          osc1(iies,iiis,iG,igain) = mean(squeeze(mean(squeeze(out.osc1),1)));
          
        end
      end
    end
  end
  
  save(sprintf('~/pmod/proc/pmod_wc_wholebrain_detosc_all_v%d.mat',vv),'osc1')
  
  for iies = 1 : length(Ies)
    iies
    for iiis = 1 : length(Iis)
      for iG = 1:length(Gg)
        for igain = 1:length(Gains)
          
          warning('Deleting...')
          delete(sprintf('~/pmod/proc/pmod_wc_wholebrain_detosc_Ie%d_Ii%d_G%d_gain%d_v%d.mat',iies,iiis,iG,igain,vv))
          delete(sprintf('~/pmod/proc/pmod_wc_wholebrain_detosc_Ie%d_Ii%d_G%d_gain%d_v%d_processing.txt',iies,iiis,iG,igain,vv))
        end
      end
    end
  end
else
  load(sprintf('~/pmod/proc/pmod_wc_wholebrain_detosc_all_v%d.mat',vv))
end

%%

vv = 6
for iies = 1: length(Ies)
  iies
  for iiis = 1: length(Iis)
    for iG = 1 : length(Gg)
      for igain = 1 : length(Gains)
%             
        try 
        load(sprintf(['~/pmod/proc/' 'pmod_wc_wholebrain_detosc_Ie%d_Ii%d_G%d_gain%d_v%d.mat'],iies,iiis,iG,igain,vv))
          catch me
          delete(sprintf(['~/pmod/proc/' 'pmod_wc_wholebrain_detosc_Ie%d_Ii%d_G%d_gain%d_v%d_processing.txt'],iies,iiis,iG,igain,vv))
          warning('!')
        end
        
      end
    end
  end
end