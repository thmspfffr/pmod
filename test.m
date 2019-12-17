%% pmod_wc_wholebrain_detosc
% Determine oscillatory regime
%-------------------------------------------------------------------------

clear
outdir = '~/pmod/proc/';

% restoredefaultpath;matlabrc

% 29-05-2018: fit E and I through resting state recordings. Then obtain
% changes in E and I due to task from recordings and keep those parameters
% fixed for the drug simulations. Vary excitability and gain for the drug
% recordings.

%-------------------------------------------------------------------------
% VERSION 1: 20-10-2018
% %-------------------------------------------------------------------------
% v           = 1;
% Ies         = -4:0.05:-1;
% Iis         = -5:0.05:-2;
% Gg          = 0:0.05:2;
% Gains       = 0; 
% nTrials     = 1;
% tmax        = 1000;  % in units of tauE
% EC          = 0;
% dt          = 0.01;
%-------------------------------------------------------------------------
% VERSION 1: 20-10-2018
% %-------------------------------------------------------------------------
v           = 111;
Ies         = -4:0.05:-1;
Iis         = -5:0.05:-2;
Gg          = 0:0.05:2;
Gains       = 0; 
nTrials     = 3;
tmax        = 1000;  % in units of tauE
EC          = 0;
dt          = 0.01;
%-------------------------------------------------------------------------

% EXCLUDE CERTAIN REGIONS - BCN ordering
k = 1 : 90;
exclude_bcn = [11 15 21 36 37 38 39 52 53 54 55 70 76 80];
include_bcn = find(~ismember(k,exclude_bcn));

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

tspan=0:dt:tmax;
L = length(tspan);

ds = 10;
Tds = length(0:ds*dt:tmax)-1;
tauEsec = 0.009; % in seconds
resol = ds*dt*tauEsec;

%%


for igain = 1
  for iG = 28
    for iies = 50
      for iiis = 49
        
        
        if ~exist(sprintf('~/pmod/proc/detosc/v%d/',v))
          mkdir(sprintf(['~/pmod/proc/detosc/' 'v%d'],v))
        end
        
        outdir = sprintf(['~/pmod/proc/detosc/' 'v%d/'],v);
        
        fn = sprintf('pmod_wc_wholebrain_detosc_Ie%d_Ii%d_G%d_gain%d_v%d',iies,iiis,iG,igain,v);
%         if tp_parallel(fn,outdir,1,0)
%           continue
%         end

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
          tic
          fprintf('Rest, Ie%d, Ii%d, trial%d ...\n',iies,iiis,tr)
          r   = 0.001*rand(2*N,1);
          R   = zeros(Tds,N);
          Ri  = zeros(Tds,N);
          tt  = 0;
          %         transient:
          for t = 1:7000
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
            [out.osc1(tr,:,i) ] = tp_detect_osc(R(:,i));    
          end
          
%           dfa = tp_dfa(R,[1 5],1/resol,0.5,20)
          toc
        end
               
%         save(sprintf([outdir '%s.mat'],fn),'out')
%         
%         while 1
%           try
%             load(sprintf([outdir '%s.mat'],fn))
%           catch me
%             save(sprintf([outdir '%s.mat'],fn),'out')
%             continue
%           end
%           break
%         end
  
%         tp_parallel(fn,outdir,0,0)

      end
    end
  end
end

error('!')

%% LOAD / DELETE ALL FILES
% 
% outdir = sprintf(['~/pmod/proc/detosc/' 'v%d/'],v);
% 
% d=dir([outdir '*.mat']);
% 
%   if numel(d)==length(Ies)*length(Iis)*length(Gg)*length(Gains)
% 
%   vv = v;
%   outdir = sprintf('~/pmod/proc/detosc/v%d/',vv);
% 
%   for iies = 1 : length(Ies)
%     iies
%     for iiis = 1 : length(Iis)
%       for iG = 1:length(Gg)
%         for igain = 1:length(Gains)
% 
%           load(sprintf([outdir 'pmod_wc_wholebrain_detosc_Ie%d_Ii%d_G%d_gain%d_v%d.mat'],iies,iiis,iG,igain,vv))
%           osc1(iies,iiis,iG,igain) = mean(squeeze(mean(squeeze(out.osc1),1)));
% 
%         end
%       end
%     end
%   end
% 
%   save(sprintf([outdir 'pmod_wc_wholebrain_detosc_all_v%d.mat'],vv),'osc1')
% 
%   % DELETE OLD FILES
%   for iies = 1 : length(Ies)
%     iies
%     for iiis = 1 : length(Iis)
%       for iG = 1:length(Gg)
%         for igain = 1:length(Gains)
% 
%           warning('Deleting...')
%           delete(sprintf([outdir 'pmod_wc_wholebrain_detosc_Ie%d_Ii%d_G%d_gain%d_v%d.mat'],iies,iiis,iG,igain,vv))
%   %         delete(sprintf([outdir 'pmod_wc_wholebrain_detosc_Ie%d_Ii%d_G%d_gain%d_v%d_processing.txt'],iies,iiis,iG,igain,vv))
%         end
%       end
%     end
%   end
% % 
% end
