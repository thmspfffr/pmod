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
[wavelet,f]=tp_mkwavelet(11.3137,0.5,(1/resol));
delta_time = 6/pi./(f(2)-f(1));
delta_time = round(delta_time*1000)/1000;
t_shift    = delta_time;
n_win      = round(delta_time*(1/resol));
n_shift    = round(t_shift*(1/resol));
nseg=floor((L/10-n_win)/n_shift+1);
% ----------------
% BAND PASS FILTER
% ----------------
flp = 9;           % lowpass frequency of filter
fhi = 13;
k     = 4;                  % 2nd order butterworth filter
fnq   = 1/(2*resol);       % Nyquist frequency
Wn    = [flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt,afilt] = butter(k,Wn);
% ----------------

%%

for igain = 6
  for iG = 6
    for iies = 43
      for iiis = 35
        
        out.FC_env    = zeros(N,N,1);
        out.COV_env    = zeros(N,N,1);
        
        if ~exist(sprintf('~/pmod/proc/numerical/v%d/',v))
          mkdir(sprintf(['~/pmod/proc/numerical/' 'v%d'],v))
        end
        
        % save configuration, with all above parameters
        if ~exist(sprintf([outdir 'pmod_wc_wholebrain_final_parameters_v%d.mat'],v))
          save(sprintf([outdir 'pmod_wc_wholebrain_final_parameters_v%d.mat'],v))
        end

        outdir = sprintf(['~/pmod/proc/numerical/v%d/'],v);
        
        fn = sprintf('pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d',iies,iiis,iG,igain,v);
%         if tp_parallel(fn,outdir,1,0)
%           continue
%         end
        
        fprintf('Gain%d, Coupling%d, Ie%d, Ii%d...\n',igain,iG,iies,iiis)
        
        tic
        g = Gg(iG);
        W = [wEE*eye(N)+g*C -wEI*eye(N); wIE*eye(N) -wII*eye(N)];
        
        % Control parameters
        %--------------------
        out.Ie   = Ies(iies);
        out.Ii   = Iis(iiis);
        out.Gain = Gains(igain);
        
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
        Io(N+1:2*N) = Iis(iiis);
        
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
          
          out.fc_FR = rc;
          
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

          out.alphapow(:,tr) = squeeze(mean(PSD(frequencies>=flp&frequencies<=fhi,:,:),1));
          
          % ----------------------------------
          % WAVELET ANALYSIS
          % ----------------------------------
          for j=1:nseg
            dloc2=rE((j-1)*n_shift+1:(j-1)*n_shift+n_win,:)';
            dataf(j,:)=abs(dloc2*wavelet).^2; 
          end
          out.rc_wl = corr(dataf);
          out.rc_wl_cov = cov(dataf);
          % ----------------------------------
          % EXTRACT PEAK FREQ
          % ---------------------------
          [~,peak_idx]=max(smooth(mean(PSD(f>4,:),2),20));
          out.peakfreq = f(peak_idx+find(f<4,1,'last'));
          clear PSD rc fc_env
         
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
for v = [3]
  for igain = 13 : length(Gains)
    igain
    for iG = 1 %: length(Gg)
      for iies = 1: length(Ies)
        iies
        for iiis = 1: length(Iis)
        
        delete( sprintf('~/pmod/proc/numerical/v%d/pmod_wc_wholebrain_final_Ie%d_Ii%d_G%d_gain%d_v%d.mat',v,iies,iiis,iG,igain,v))

        end
%           fn = sprintf('/home/tpfeffer/pmod/proc/pmod_wc_wholebrain_detosc_Ie%d_Ii%d_G%d_gain%d_v%d.mat',iies,iiis,iG,igain,v);
%           delete(fn)
%         end
      end
    end
  end
end

% %%
% mask    = logical(triu(ones(76,76),1));
% for iies= 1 : 31
%   for iiis = 1 : 31
%     
%     load(sprintf('~/pmod/proc/numerical/v3/pmod_wc_wholebrain_final_Ie%d_Ii%d_G1_gain1_v3.mat',iies,iiis))
%     aaa(iies,iiis) = mean(out.FC_env(mask));
%   end
% end


