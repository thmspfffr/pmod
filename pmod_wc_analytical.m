%% pmod_wc_analytical

% Computes correlations between two Wilson-Cowan nodes analytically
% Code from Adrian Ponce-Alvarez, edits Tom Pfeffer


clear all;
outdir  = '~/pmod/proc/';
% -------------------------------
% VERSION 1 
% -------------------------------
v = 1;
Ies = -5:.1:3; %-4:.1:-1;
Iis = -13:.1:2; %-4:.1:-1;
all_gains = 0:0.1:0.6;
% -------------------------------
% v = 3;
% Ies = -4:0.025:-1;
% Iis = -5:0.025:-2;
% all_gains = 0:0.1:0.6;
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

% Connectivity:
W11 = [wEE -wEI; wIE -wII];
W22 = W11;
W12 = [g 0; 0 0];
W21 = [g 0; 0 0];

W = [W11 W12; W21 W22];

p = 1.5;

numIes = length(Ies);
numIis = length(Iis);
%%
for igain = 1 : length(all_gains)
  
  Gain = all_gains(igain);
  if ~exist(sprintf([outdir 'pmod_wc_analytical_gain%d_v%d_processing.txt'],igain,v))
    system(['touch ' outdir sprintf('pmod_wc_analytical_gain%d_v%d_processing.txt',igain,v)]);
  else
    continue
  end
  
  fprintf('Gain %d...\n',igain)
  for Q = 2  % conditions (no gain change, equivalent gain change, daE<daI, daE>daI)
    
    % transfer functions:
    % gains are given by 1/aE and 1/aI
    if Q == 1
      gE = 1;
      gI = 1;
    elseif Q == 2
      gE = 1+Gain;
      gI = 1+Gain;
    elseif Q == 3
      gE = 1+p*Gain;
      gI = 1+Gain;
    elseif Q == 4
      gE = 1+Gain;
      gI = 1+p*Gain;
    end
    
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
          
%           display([Q Ie Ii maxR-minR])
          
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
          
%           display(Corr)
          out.CeE(k,l) = Corr(1,3);
          
        else
          
          out.CeE(k,l) = NaN;
          
        end
        
        
      end
    end
  end
  
  save(sprintf([outdir 'pmod_wc_analytical_gain%d_v%d.mat'],igain,v),'out')

end
% save WCs_Correlations_ParamSpace CeE p Gain W g tau Iis Ies
%%
v = 1;
for igain = 1 : 7
  
  load(sprintf([outdir 'pmod_wc_analytical_gain%d_v%d.mat'],igain,v))

  CeE(:,:,igain) = out.CeE;
  
end

%%
y=[];
x=[];
for j=1:size(CeE,2);
  temp=find(isnan(CeE(:,j,1)),1,'first');
  if temp<size(CeE,2);
    y(j)=Ies(temp);
    x(j)=Iis(j);
  end
end

load redblue.mat
idxE = find(Ies==-4,1,'first'):find(Ies==-1,1,'first');
idxI = find(Iis==-5,1,'first'):find(Iis==-2,1,'first');
% CeE(isnan(CeE)) = 0;
figure; set(gcf,'color','w');

a = subplot(3,2,1); hold on
imagescnan(Iis,Ies,CeE(:,:,1),[0 .25]); colormap(jet)
xlabel('Background input (to I)'); ylabel('Background input (to E)');
tp_editplots; colormap(a,'plasma')
axis square tight

labI = Iis(idxI);
labE = Ies(idxE);

b = subplot(3,2,2); hold on

imagescnan(CeE(idxE,idxI,1),[0 .25])
plot(find(Iis==-3.7)-idxI(1),find(round(Ies*10)==-2.7*10)-idxE(1),'ko','markersize',5,'markerfacecolor','r')
    plot(find(round(Iis*10)==-3.1*10)-idxI(1),find(round(Ies*10)==-2.1*10)-idxE(1),'ko','markersize',5,'markerfacecolor','y')
  
xlabel('Background input (to I)'); ylabel('Background input (to E)');
tp_editplots; colormap(b,plasma)
axis tight equal
set(gca,'YTick',1:10:length(labE),'YTickLabels',num2cell(labE(1:10:end)))
set(gca,'XTick',1:10:length(labI),'XTickLabels',num2cell(labI(1:10:end)))

if v==1
  idxE = find(Ies==-4,1,'first'):find(Ies==-1,1,'first');
  idxI = find(Iis==-5,1,'first'):find(Iis==-2,1,'first');
  
  for igain = [2 2 4 7]
    if igain== 7
      c = subplot(3,2,6); hold on 
    else
      c = subplot(3,2,1+igain); hold on 
    end
    imagescnan(CeE(idxE,idxI,igain)-CeE(idxE,idxI,1),[-.1 .1])
%     if igain == 
    plot(find(Iis==-3.7)-idxI(1),find(round(Ies*10)==-2.7*10)-idxE(1),'ko','markersize',5,'markerfacecolor','r')
    plot(find(round(Iis*10)==-3.1*10)-idxI(1),find(round(Ies*10)==-2.1*10)-idxE(1),'ko','markersize',5,'markerfacecolor','y')
    xlabel('Background input (to I)'); ylabel('Background input (to E)');
    tp_editplots; colormap(c,redblue)
    axis tight equal square
    labI = Iis(idxI);
    labE = Ies(idxE);
    set(gca,'YTick',1:10:length(labE),'YTickLabels',num2cell(labE(1:10:end)))
    set(gca,'XTick',1:10:length(labI),'XTickLabels',num2cell(labI(1:10:end)))

    set(gcf,'renderer','Painters')
  end
end
print(gcf,'-depsc2',sprintf('~/pmod/plots/pmod_wc_analytical_v%d.eps',v))
%%
igain = 3;

figure; set(gcf,'color','w');
c = subplot(3,2,1); hold on 

imagescnan(CeE(idxE,idxI,igain)-CeE(idxE,idxI,1),[-.1 .1])
plot(find(Iis==-3.7)-idxI(1),find(round(Ies*10)==-2.7*10)-idxE(1),'ko','markersize',5,'markerfacecolor','r')
plot(find(round(Iis*10)==-3.1*10)-idxI(1),find(round(Ies*10)==-2.1*10)-idxE(1),'ko','markersize',5,'markerfacecolor','y')
xlabel('Background input (to I)'); ylabel('Background input (to E)');
tp_editplots; colormap(c,redblue)
axis tight square
labI = Iis(idxI);
labE = Ies(idxE);
set(gca,'YTick',1:10:length(labE),'YTickLabels',num2cell(labE(1:10:end)))
set(gca,'XTick',1:10:length(labI),'XTickLabels',num2cell(labI(1:10:end)))

c = subplot(3,2,2); hold on 

imagescnan((CeE(idxE,idxI,igain)-CeE(idxE,idxI,1))./CeE(idxE,idxI,1),[-.85 .85])
plot(find(round(Iis*10)==-4.3*10)-idxI(1),find(round(Ies*10)==-3.3*10)-idxE(1),'ko','markersize',5,'markerfacecolor','r')
plot(find(round(Iis*10)==-3.7*10)-idxI(1),find(round(Ies*10)==-2.7*10)-idxE(1),'ko','markersize',5,'markerfacecolor','y')
% xlabel('Background input (to I)'); ylabel('Background input (to E)');
tp_editplots; colormap(c,redblue)
axis tight square
labI = Iis(idxI);
labE = Ies(idxE);
set(gca,'YTick',1:10:length(labE),'YTickLabels',num2cell(labE(1:10:end)))
set(gca,'XTick',1:10:length(labI),'XTickLabels',num2cell(labI(1:10:end)))


r = 100*(CeE(idxE,idxI,igain)-CeE(idxE,idxI,1))./(CeE(idxE,idxI,1));

idxI_rest = find(round(Iis*10)==-3.7*10)-idxI(1);
idxE_rest = find(round(Ies*10)==-2.7*10)-idxE(1);
idxI_task = find(round(Iis*10)==-3.1*10)-idxI(1);
idxE_task = find(round(Ies*10)==-2.1*10)-idxE(1);

subplot(2,2,3); hold on
bar([1 2 3],[r(idxE_rest,idxI_rest) r(idxE_task,idxI_task) r(idxE_rest,idxI_rest)-r(idxE_task,idxI_task)]); axis square
axis([0.5 3.5 -60 20]); tp_editplots; ylabel('Change in correlation [in %]')

idxI_rest = find(round(Iis*10)==-4.3*10)-idxI(1);
idxE_rest = find(round(Ies*10)==-3.3*10)-idxE(1);
idxI_task = find(round(Iis*10)==-3.7*10)-idxI(1);
idxE_task = find(round(Ies*10)==-2.7*10)-idxE(1);

subplot(2,2,4); hold on
bar([1 2 3],[r(idxE_rest,idxI_rest) r(idxE_task,idxI_task) r(idxE_rest,idxI_rest)-r(idxE_task,idxI_task) ]); axis square
axis(gca,[0.5 3.5 -60 20]); tp_editplots; ylabel('Change in correlation [in %]')

print(gcf,'-depsc2',sprintf('~/pmod/plots/pmod_wc_analytical_bar_v%d.eps',v))

%%

fc1 = rand(50,20,2,100);
fc2 = rand(50,20,2,100);

scatter(squeeze(mean(mean(fc1(:,:,1,:)))),squeeze(mean(mean(fc1(:,:,2,:))))-squeeze(mean(mean(fc1(:,:,1,:)))))

h = ttest(fc1(:,:,1,:),fc1(:,:,2,:),'dim',2);
sum(squeeze(h))

