%% pmod_wc_analytical

% Computes correlations between two Wilson-Cowan nodes analytically
% Code from Adrian Ponce-Alvarez, edits Tom Pfeffer

clear all;

% outdir  = '~/pmod/proc/analytical/';

% -------------------------------
% VERSION 1
% -------------------------------
v = 3;
Ies = -4:0.025:-1;
Iis = -5:0.025:-2;
all_gains = 0:0.02:0.3;
all_inh = 1:-0.05:0.7;

% -------------------------------

% fixed params:
N = 4;

% Connectivity:
wII=4;
wIE=16;
wEI=12;
wEE=12;
g = 1;

tauE = 1;
tauI = 2; %1
tau = [tauE;tauI;tauE;tauI];

Io=zeros(N,1);

dt=0.001;

tmax = 20; %5000;
tspan=0:dt:tmax;
L = length(tspan);
ds = 5;
Tds = length(0:ds*dt:tmax)-1;
sigma = 0.0005;
Qn = (sigma*dt)^2*eye(N);

p = 1.5;

numIes = length(Ies);
numIis = length(Iis);

%%

for igain = 1 : length(all_gains)
  
  Gain = all_gains(igain);
  
  if ~exist(sprintf(['~/pmod/proc/analytical/v%d'],v))
    mkdir(sprintf(['~/pmod/proc/analytical/v%d'],v))
  end
  outdir = sprintf(['~/pmod/proc/analytical/v%d/'],v);
  
  if ~exist(sprintf([outdir 'pmod_wc_analytical_gain%d_v%d_processing.txt'],igain,v))
    system(['touch ' outdir sprintf('pmod_wc_analytical_gain%d_v%d_processing.txt',igain,v)]);
  else
    continue
  end
  
  for iinh = 1 : length(all_inh)
    % Connectivity:
    W11 = [wEE -wEI*all_inh(iinh); wIE -wII];
    W22 = W11;
    W12 = [g 0; 0 0];
    W21 = [g 0; 0 0];

    W = [W11 W12; W21 W22];


    fprintf('Gain %d...\n',igain)

    gE = 1+Gain;
    gI = 1+Gain;
    aE = 1/gE;
    aI = 1/gI;

    Fe = @(x) 1./(1 + exp(-x/aE) );
    Fi = @(x) 1./(1 + exp(-x/aI) );
    F = @(x) [feval(Fe,x(1));feval(Fi,x(2));feval(Fe,x(3));feval(Fi,x(4))];

    for l = 1:numIis

      l

      Ii = Iis(l);
      Io(2) = Ii;
      Io(4) = Ii;

      for k = 1:numIes

        Ie = Ies(k);
        Io(1) = Ie;
        Io(3) = Ie;

        r = rand(N,1);
        R = zeros(Tds,N);
        tt = 0;

        % get fixed points:
        for t = 1:300000
          u = W*r + Io;
          K = feval(F,u);
          r = r + dt*(-r + K)./tau;
        end

        for t = 1:L
          u = W*r + Io;
          K = feval(F,u);
          r = r + dt*(-r+ K)./tau;
          if mod(t,ds)==0
            tt=tt+1;
            R(tt,:)=r;
          end
        end

        RstatE = R(:,1);
        maxR = max(RstatE);
        minR = min(RstatE);

        if maxR-minR<10e-10

          Cov = zeros(N);

          for t = 1:10000

            u = W*r + Io;
            K = feval(F,u);
            r = r + dt*(-r+ K)./tau;

            re = r(1);
            ri = r(2);
            rE = r(3);
            rI = r(4);

            % Jacobian:
            Aee = -1 + wEE*1/aE*re*(1-re);
            Aei = -wEI*1/aE*re*(1-re);
            Aie = wIE*1/aI*ri*(1-ri)*(tauE/tauI);
            Aii = -1*(tauE/tauI)- wII*1/aI*ri*(1-ri)*(tauE/tauI);

            AeE = g*1/aE*re*(1-re);
            AeI = 0;
            AiE = 0;
            AiI = 0;

            AEe = g*1/aE*rE*(1-rE);
            AIe = 0;
            AEi = 0;
            AIi = 0;

            AEE = -1/tauE + wEE*1/aE*rE*(1-rE);
            AEI = -wEI*1/aE*rE*(1-rE);
            AIE = wIE*1/aI*rI*(1-rI)*(tauE/tauI);
            AII = -1*(tauE/tauI) - wII*1/aI*rI*(1-rI)*(tauE/tauI);

            Jmat = [Aee Aei AeE AeI; ...
              Aie Aii AiE AiI; ...
              AEe AEi AEE AEI; ...
              AIe AIi AIE AII];

            delta = (Jmat*Cov + Cov*Jmat' + Qn)*dt;
            Cov = Cov + delta;

            tol = norm(delta);
            if tol<10e-18
              break
            end
          end

          Corr=zeros(N);
          for m=1:N
            for n=1:N
              Corr(m,n)=Cov(m,n)/( sqrt(Cov(m,m))*sqrt(Cov(n,n)) );
            end
          end

          out.CeE(k,l,iinh) = Corr(1,3);

        else
          out.CeE(k,l,iinh) = NaN;
        end
      end
    end
  
  
    save(sprintf([outdir 'pmod_wc_analytical_gain%d_v%d.mat'],igain,v),'out')
  end
end

% save WCs_Correlations_ParamSpace CeE p Gain W g tau Iis Ies

%%
Ies = -4:0.025:-1;
Iis = -5:0.025:-2;

v = 3;
outdir = sprintf('~/pmod/proc/analytical/v%d/',v);

for igain = 1:11
  load(sprintf([outdir 'pmod_wc_analytical_gain%d_v%d.mat'],igain,v))
  CeE(:,:,igain) = out.CeE(:,:,1);
end

%%

load redblue.mat

for igain = 11
  
  h=figure; set(gcf,'color','w');
  
  b=subplot(2,3,1); hold on
  %   igain = gain;
  par = CeE(:,:,igain)-CeE(:,:,1);
  
  c=imagesc(par,[-.05 .05]); axis square
  set(c,'AlphaData',~isnan(par))
  
  
  colormap(b,redblue)
  axis square; set(gca,'ydir','normal')
  axis([1 121 1 121])
  set(gca,'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  set(gca,'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
  
  xlabel('Background input to I'); ylabel('Background input to E'); tp_editplots
  hf=text(10,60,sprintf('Oscillations'),'fontsize',8,'fontangle','italic','color',[.5 .5 .5],'parent',b)
  set(hf,'Rotation',45)
  
  ii_rest_pbo = -3.7; ii_task_pbo = -3.3;
  ie_rest_pbo = -2.7; ie_task_pbo = -2.3;
  ii_rest_dpz = -4.1; ii_task_dpz = -3.7;
  ie_rest_dpz = -3.1; ie_task_dpz = -2.7;

  scatter(find(round(Iis*1000)==round(ii_rest_dpz*1000)),find(round(Ies*1000)==round(ie_rest_dpz*1000)),15,'o','markeredgecolor','k','markerfacecolor',[.8 .8 .8])
  scatter(find(round(Iis*1000)==round(ii_task_dpz*1000)),find(round(Ies*1000)==round(ie_task_dpz*1000)),15,'o','markeredgecolor','k','markerfacecolor',[1 1 1])
  
  corr_rest_pbo=CeE(find(round(Ies.*1000)==ie_rest_pbo*1000),find(round(Iis*1000)==ii_rest_pbo*1000),1)
  corr_rest_dpz=CeE(find(round(Ies.*1000)==ie_rest_dpz*1000),find(round(Iis*1000)==ii_rest_dpz*1000),igain)
  corr_task_pbo=CeE(find(round(Ies.*1000)==ie_task_pbo*1000),find(round(Iis*1000)==ii_task_pbo*1000),1)
  corr_task_dpz=CeE(find(round(Ies.*1000)==ie_task_dpz*1000),find(round(Iis*1000)==ii_task_dpz*1000),igain)
  
  subplot(2,3,2); hold on
  
  bar(1,corr_rest_pbo,'facecolor',[.8 .8 .8])
  bar(2,corr_rest_dpz,'facecolor',[1 1 1])
  bar(3,corr_task_pbo,'facecolor',[1 0 0])
  bar(4,corr_task_dpz,'facecolor',[1 .5 .5])
  
  axis square; tp_editplots
  axis([0 5 0 0.3]);
  set(gca,'XTick',[1 2 3 4],'XTickLabels',{'Rest (Pbo)';'Rest (Atx)';'Task (Pbo)';'Task (Atx)'})
  xtickangle(gca,45)
  
  b=subplot(2,3,4); hold on
  par = CeE(:,:,igain)-CeE(:,:,1);
  c=imagesc(par,[-.05 .05]); axis square
  set(c,'AlphaData',~isnan(par))
  colormap(b,redblue);
  
%   ii_rest = -3.7; ii_task = -3.3;
%   ie_rest = -2.7; ie_task = -2.3;
  
  scatter(find(round(Iis*1000)==round(ii_rest_pbo*1000)),find(round(Ies*1000)==round(ie_rest_pbo*1000)),15,'o','markeredgecolor','k','markerfacecolor',[.8 .8 .8])
  scatter(find(round(Iis*1000)==round(ii_task_pbo*1000)),find(round(Ies*1000)==round(ie_task_pbo*1000)),15,'o','markeredgecolor','k','markerfacecolor',[1 1 1])
  
%   colormap(b,plasma); hold on
  axis square; set(gca,'ydir','normal')
  set(gca,'XTick',1:40:length(Iis),'XTickLabels',num2cell(Iis(1:40:end)))
  set(gca,'YTick',1:40:length(Ies ),'YTickLabels',num2cell(Ies(1:40:end)))
  xlabel({'Background input to I'}); ylabel({'Background input to E'}); tp_editplots
  axis([1 121 1 121])
  
  corr_rest_pbo=CeE(find(Ies==ie_rest_pbo),find(Iis==ii_rest_pbo),1)
  corr_rest_atx=CeE(find(Ies==ie_rest_pbo),find(Iis==ii_rest_pbo),igain)
  corr_task_pbo=CeE(find(round(Ies.*1000)==ie_task_pbo*1000),find(round(Iis.*1000)==ii_task_pbo*1000),1)
  corr_task_atx=CeE(find(round(Ies.*1000)==ie_task_pbo*1000),find(round(Iis.*1000)==ii_task_pbo*1000),igain)
  
  
  subplot(2,3,5); hold on
  
  bar(1,corr_rest_pbo,'facecolor',[.8 .8 .8])
  bar(2,corr_rest_atx,'facecolor',[1 1 1])
  bar(3,corr_task_pbo,'facecolor',[1 0 0])
  bar(4,corr_task_atx,'facecolor',[1 .5 .5])
  
  axis square; tp_editplots
  axis([0 5 0 0.3]);
  set(gca,'XTick',[1 2 3 4],'XTickLabels',{'Rest (Pbo)';'Rest (Atx)';'Task (Pbo)';'Task (Atx)'})
  xtickangle(gca,45)
  
  
  
end

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_wc_analytical_gain%d_v%d.pdf',igain,v))

%% TEST ALTERNATIVE MECHANISMS SYSTEMATICALLY
