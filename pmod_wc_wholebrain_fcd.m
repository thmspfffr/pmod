%% pmod_wc_wholebrain_fcd
% Stochastic simulation of 2*N WC nodes during "rest"
%-------------------------------------------------------------------------- 

clear

%--------------------------------------------------------------------------
% VERSION 1: DETERMINE GLOBAL COUPLING PARAMETER
% %------------------------------------------------------------------------
v           = 1;
Ies         = -4:0.05:-1;
Iis         = -5:0.05:-2;
Gg          = 0.1:0.1:1.5;
Gains       = 0;
nTrials     = 1;
tmax        = 67000;  % in units of tauE
EC          = 0;
% -------------------------------------------------------------------------
% VERSION 2: SAME AS ABOVE, BUT INCLUDING FCD ON FIRING RATES
% %------------------------------------------------------------------------
v           = 2;
Ies         = -4:0.05:-1;
Iis         = -5:0.05:-2;
Gg          = 0.1:0.1:1.5;
Gains       = 0;
nTrials     = 1;
tmax        = 67000;  % in units of tauE
EC          = 0;
% -------------------------------------------------------------------------

% EXCLUDE CERTAIN REGIONS - BCN layout
% ------------------
k = 1 : 90;
exclude_bcn = [11 15 21 36 37 38 39 52 53 54 55 70 76 80];
include_bcn = find(~ismember(k,exclude_bcn));
% ------------------

% load connectome
load ~/sc90.mat %Bea SC
C = SC;
C = C/max(C(C>0));
C = C(include_bcn,include_bcn);
N = size(C,1);

addpath ~/pconn/matlab
%--------------------------------------------------------------------------
% PARAMETER DEFINITIONS
%--------------------------------------------------------------------------

% Connectivity:
wII   = 4;
wIE   = 16;
wEI   = 12;
wEE   = 12;
% time constants
tauE          = 1;
tauI          = 2;
tau           = zeros(2*N,1);
tau(1:N)      = tauE;
tau(N+1:2*N)  = tauI;

dt=0.01;
tspan=0:dt:tmax;
L = length(tspan);
clear tspan

ds        = 10;
Tds       = length(0:ds*dt:tmax)-1;
tauEsec   = 0.009; % in seconds
resol     = ds*dt*tauEsec;

sigma = 0.0005;

isub = find( triu(ones(N)) - eye(N) );
mask      = logical(tril(ones(76,76),-1));

% ----------------
% WAVELETS
% ----------------
[wavelet,f]=tp_mkwavelet(11.3137,0.5,(1/resol));
delta_time = 6/pi./(f(2)-f(1));
delta_time = round(delta_time*1000)/1000;
t_shift    = delta_time/2;
n_win      = round(delta_time*(1/resol));
n_shift    = round(t_shift*(1/resol));
nseg=floor((L/10-n_win)/n_shift+1);

%%
for igain = 1 : length(Gains)
  for iG = 1 : length(Gg)
    for iies = 1: length(Ies)
      for iiis = 1: length(Iis)

        if ~exist(sprintf('~/pmod/proc/fcd/v%d/',v))
          mkdir(sprintf(['~/pmod/proc/fcd/' 'v%d'],v));
        end
        outdir = sprintf(['~/pmod/proc/fcd/v%d/'],v);

        fn = sprintf('pmod_wc_wholebrain_fcd_Ie%d_Ii%d_G%d_gain%d_v%d',iies,iiis,iG,igain,v);
        if tp_parallel(fn,outdir,1,0)
          continue
        end
        
        fprintf('Gain%d, Coupling%d, Ie%d, Ii%d...\n',igain,iG,iies,iiis)
        
        tic
        g = Gg(iG);
        W = [wEE*eye(N)+g*C -wEI*eye(N); wIE*eye(N) -wII*eye(N)];
        %--------------------
        % Control params.
        %--------------------
        out.Ie = Ies(iies);
        out.Ii = Iis(iiis);
        out.Gain = Gains(igain);
        out.coupl = g;
        
        % Working point:
        Io          = zeros(2*N,1);
        Io(1:N)     = out.Ie;
        Io(N+1:2*N) = out.Ii;
        
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
        PSD         = zeros(length(pp),N,nTrials);
        frequencies = freq100(1:10:end)';
        
        % RUN SIMULATION
        % ---------------------
        for tr=1:nTrials
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
          
          % FC matrix
          % ---------------------
          rc       	= corrcoef(rE);
          fc      	= rc(isub);
          
          for i=1:N
            
            % COMPUTE POWER SPECTRUM
            % ---------------------------
            f           = rE(:,i) - mean(rE(:,i));
            xdft        = fft(f);
            xdft        = xdft(1:floor(Tds/2)+1);
            pw          = (1/(Tds/2)) * abs(xdft).^2;
            psd         = pw(freqs<100 & freqs>1);
            f           = freqs(freqs<100 & freqs>1);
            fnew        = f(1:10:end);
            psd         = psd(1:10:end);
            PSD(:,i,tr) = psd';
            f           = fnew;
            
          end
          
          [wavelet,~,opt] = tp_mkwavelet(11.3137,0.5,(1/resol));
          nseg=floor((L/10-opt.n_win)/opt.n_shift+1);
%         
%           out.alphapow(:,tr) = squeeze(mean(PSD(frequencies>=flp&frequencies<=fhi,:,:),1));
          
          % EXTRACT PEAK FREQ
          % ---------------------------
          [~,peak_idx]=max(smooth(mean(PSD(f>4,:),2),20));
          out.peakfreq = f(peak_idx+find(f<4,1,'last'));
                    
          for j=1:nseg
            dloc2=rE((j-1)*n_shift+1:(j-1)*n_shift+n_win,:)';
            dataf(j,:)=abs(dloc2*wavelet).^2;
            seg_data(:,:,j) = dloc2;
          end
          out.rc_wl = corr(dataf);
          out.rc_wl_cov = cov(dataf);
          
          % 60s window length, 40s overlap (Deco et al., 2017, Sci Rep)
          segleng = floor(40/(size(wavelet,1)*(resol)));
          segshift = floor(20/(size(wavelet,1)*(resol)));
          nseg = floor((size(dataf,1)-segleng)/segshift+1);
                    
          for iseg = 1 : nseg
            dloc = dataf((iseg-1)*segshift+1:(iseg-1)*segshift+segleng,:);
            r_env(:,:,iseg) = corrcoef(dloc);
            dloc = seg_data(:,:,(iseg-1)*segshift+1:(iseg-1)*segshift+segleng);
            dloc = reshape(dloc,[76 opt.n_win*segleng]);
            r_env_fr(:,:,iseg)=corrcoef(dloc');
          end
          
          
          for iseg = 1 : nseg
            for jseg = 1 : nseg
              tmp1 = r_env(:,:,iseg);
              tmp2 = r_env(:,:,jseg);
              out.fcd_env(iseg,jseg) = corr(tmp1(mask),tmp2(mask));
              tmp1 = r_env_fr(:,:,iseg);
              tmp2 = r_env_fr(:,:,jseg);
              out.fcd_env_fr(iseg,jseg) = corr(tmp1(mask),tmp2(mask));
            end
          end
          
          % COMPUTE FCD BASED ON FIRING RATES
          % ---------------------------
          

          
          toc
        end
        
        save([outdir sprintf('%s.mat',fn)],'out')
        
%         make sure file is saved
        while 1
          try
            load([outdir sprintf('%s.mat',fn)])
            break
          catch me
            save([outdir sprintf('%s.mat',fn)],'out')
          end
        end
        
        clear out dataf nseg r_env
        
        tp_parallel(fn,outdir,0)
        
      end
    end
  end
end

error('!')
%%
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

figure; set(gcf,'color','w'); hold on

subplot(1,2,1);
imagesc(nanmean(fc(:,:,:,1,2,6),3),[0.01 0.1]); axis square off
colormap(plasma)

subplot(1,2,2);
imagesc(nanmean(fc(:,:,:,2,2,6),3),[0.01 0.1]); axis square off
colormap(plasma)
%
% subplot(1,2,3);
% imagesc(nanmean(fc(:,:,:,3,2,6),3),[0.01 0.1]); axis square off
% colormap(plasma)

