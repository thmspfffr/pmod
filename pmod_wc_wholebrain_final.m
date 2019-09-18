%% pmod_wc_wholebrain_final
% Stochastic simulation of 2*N WC nodes during "rest"
%--------------------------------------------------------------------------

clear
restoredefaultpath;matlabrc

outdir = '~/pmod/proc/';

% 29-05-2018: fit E and I through resting state recordings. Then obtain
% changes in E and I due to task from recordings and keep those parameters
% fixed for the drug simulations. Vary excitability and gain for the drug
% recordings.
%--------------------------------------------------------------------------
% VERSION 1: 20-10-2018: DETERMINE GLOBAL COUPLING PARAMETER
% %------------------------------------------------------------------------
% v           = 1;
% Ies         = -4:0.1:-1;
% Iis         = -5:0.1:-2;
% Gg          = 0:0.15:3;
% Gains       = 0;
% nTrials     = 1;
% tmax        = 6500;  % in units of tauE
% EC          = 0;
%--------------------------------------------------------------------------
% VERSION 2: 20-10-2018: DETERMINE GLOBAL COUPLING PARAMETER
% %------------------------------------------------------------------------
% v           = 2;
% Ies         = -4:0.025:-1;
% Iis         = -5:0.025:-2;
% Gg          = 1.65;
% Gains       = [0 0.025:0.025:0.4 -0.025:-0.025:-0.4 0.425:0.025:0.6 0.625:0.025:0.7];
% nTrials     = 1;
% tmax        = 6500;  % in units of tauE
%--------------------------------------------------------------------------
% VERSION 2: 20-10-2018: DETERMINE GLOBAL COUPLING PARAMETER
% %------------------------------------------------------------------------

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
dt=0.01;

tspan=0:dt:tmax;
L = length(tspan);
clear tspan

ds = 10;
Tds = length(0:ds*dt:tmax)-1;
tauEsec = 0.009; % in seconds
resol = ds*dt*tauEsec;

isub = find( triu(ones(N)) - eye(N) );

%%

for igain = 1 : length(Gains)
  for iG = 1 : length(Gg)
    for iies = 1: length(Ies)
      iies
      
      if ~exist(sprintf('~/pmod/proc/numerical/v%d/',v))
        mkdir(sprintf(['~/pmod/proc/numerical/' 'v%d'],v))
      end
      
      outdir = sprintf(['~/pmod/proc/numerical/v%d/'],v);

      fn = sprintf('pmod_wc_wholebrain_final_Ie%d_G%d_gain%d_v%d',iies,iG,igain,v);
      if tp_parallel(fn,outdir,1,0)
        continue
      end
      
      % save configuration, with all above parameters
      if ~exist(sprintf([outdir 'pmod_wc_wholebrain_final_parameters_v%d.mat'],v))
        save(sprintf([outdir 'pmod_wc_wholebrain_final_parameters_v%d.mat'],v))
      end
      
      for iiis = 1: length(Iis)
        
        fprintf('Gain%d, Coupling%d, Ie%d, Ii%d...\n',igain,iG,iies,iiis)

        tic
        g = Gg(iG);
        W = [wEE*eye(N)+g*C -wEI*eye(N); wIE*eye(N) -wII*eye(N)];
        
        out(iiis).Cee       = zeros(1,1);
        out(iiis).CeeSD     = zeros(2,1);
        out(iiis).FC_env    = zeros(N,N,1);
        out(iiis).Cee_env   = zeros(1,1);
        out(iiis).CeeSD_env = zeros(2,1);
        
        out(iiis).KOPsd   = zeros(nTrials,1);
        out(iiis).KOPmean = zeros(nTrials,1);
        
        % Control params.
        %--------------------
        out(iiis).Ie = Ies(iies);
        out(iiis).Ii = Iis(iiis);
        out(iiis).Gain = Gains(igain);
        
        % Working point:
        Io=zeros(2*N,1);
        Io(1:N) = out(iiis).Ie;
        Io(N+1:2*N) = out(iiis).Ii;
        
        % transfer function:
        gE = 1 + Gains(igain);
        gI = 1 + Gains(igain);
        aE  = 1/gE;
        aI  = 1/gI;
        Fe  = @(x) 1./(1 + exp(-x/aE) );
        Fi  = @(x) 1./(1 + exp(-x/aI) );
        F   = @(x) [feval(Fe,x(1:N));feval(Fi,x(N+1:2*N))];
        
        T           = Tds*resol; %% define time of interval
        freqs       = (0:Tds/2)/T; %% find the corresponding frequency in Hz
        freq100     = freqs(freqs<100 & freqs>1);
        pp          = 1:10:length(freq100);
        PSD     = zeros(length(pp),N,nTrials);
        frequencies = freq100(1:10:end)';
        
        % RUN SIMULATION
        % ---------------------
        
        for tr=1:nTrials
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
          
%           clear ku KOP z
          
          % FC matrix
          % ---------------------
          rc       	= corrcoef(rE);
          %           out.FC  	= single(out.FC) + single(rc/nTrials);
          fc      	= rc(isub);
          flp = 9;           % lowpass frequency of filter
          fhi = 13;
          
          for i=1:N
            
            %           % COMPUTE POWER SPECTRUM
            %           % ---------------------------
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
          
          out(iiis).alphapow(:,tr) = squeeze(mean(PSD(frequencies>=flp&frequencies<=fhi,:,:),1));
          
          % EXTRACT PEAK FREQ
          % ---------------------------
          [~,peak_idx]=max(smooth(mean(PSD(f>3,:),2),20));
          out(iiis).peakfreq = f(peak_idx+find(f<4,1,'last'));
          clear PSD
          
          delt  = resol;            % sampling interval
          k     = 4;                  % 2nd order butterworth filter
          fnq   = 1/(2*delt);       % Nyquist frequency
          Wn    = [flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
          [bfilt,afilt] = butter(k,Wn);
          
          % COMPUTE DFA, EXTRACT HURST
          % ---------------------------
          env = abs(hilbert(filtfilt(bfilt,afilt,rE)));
          
          % COMPUTE CORRELATIONS BASED ON ENV
          % ---------------------------
          % Alpha range
          % ---------------------------
          rc                    = corrcoef(env.^2);
          out(iiis).FC_env      = single(out(iiis).FC_env) + single(rc/nTrials);
          fc_env                = rc(isub);
          out(iiis).Cee_env     = out(iiis).Cee_env + mean(fc)/nTrials;
          out(iiis).CeeSD_env   = out(iiis).CeeSD_env + var(fc)/nTrials;
          % ---------------------------
          clear rc fc_env

          clear rE rI env
          toc
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
  end
end

% error('!')
%%

% fprintf('WAAAAAT')
%
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
%  /home/tpfeffer/pmod/proc/pmod_wc_wholebrain_final_Ie14_Ii56_G1_gain43_v2.mat
%% DELET EFILES
% for v = [1 2 3 22]
%   for igain = 1 : length(Gains)
%     igain
%     for iG = 1 : length(Gg)
%       for iies = 1: length(Ies)
% %         pmod_wc_wholebrain_final_Ie14_Ii56_G1_gain36_v
%         iies
%         delete( sprintf('/home/tpfeffer/pmod/proc/pmod_wc_wholebrain_detosc_Ie%d_Ii*_G%d_gain%d_v%d.mat',iies,iG,igain,v))
% 
% %         for iiis = 1: length(Iis)
% %           fn = sprintf('/home/tpfeffer/pmod/proc/pmod_wc_wholebrain_detosc_Ie%d_Ii%d_G%d_gain%d_v%d.mat',iies,iiis,iG,igain,v);
% %           delete(fn)
% %         end
%       end
%     end
%   end
% end
