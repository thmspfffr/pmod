%% BIFURCATION ANALYSIS
% Computes bifurcation diagrams

clear

% fixed params:
v = 1;
% Ies = -4:0.025:-1;
Ies = -1.4;

Iis = -5:0.025:-2;
% Ies         = -4:0.025:-1;
% Iis         = -5:0.025:-2;
N = 2;
% Connectivity:
wII=4;
wIE=16;
wEI=12;
wEE=12;

tauE = 1;
tauI = 2;
% tau = zeros(2*N,1);
tau(1) = tauE;
tau(2) = tauI;
tau = tau';

Io=zeros(N,1);

tmax	= 1000; % in units of tauE
dt=0.01;
tspan=0:dt:tmax;
L = length(tspan);
clear tspan

ds = 10;
Tds = length(0:ds*dt:tmax)-1;

% transfer functions:
aE = 1;
Fe = @(x) 1./(1 + exp(-x/aE) );

aI = 1;
Fi = @(x) 1./(1 + exp(-x/aI) );

F = @(x) [feval(Fe,x(1));feval(Fi,x(2))];

% inverse functions:
FinvE = @(x) log( aE*x./(1-aE*x) );
FinvI = @(x) log( aI*x./(1-aI*x) );

numIes = length(Ies);

Re  = zeros(numIes,2);
Tr  = zeros(numIes,1);
De  = zeros(numIes,1);
Eig = zeros(numIes,N);
%%
for k = 1:numIes
  for j = 1 : length(Iis)
    
    fprintf('#%d ..\n',k)
    
    Io(1) = Ies(k);
    Io(2) = Iis(j);
    
    W = [wEE -wEI; wIE -wII];
        
    % fixed points (i.e., dx/dt = 0)
    S = @(x) -x./tau + feval(F,W*x+Io);
    x0 = .5*ones(N,1);
    fp = fsolve(S,x0);
    re = fp(1);
    ri = fp(2);
    
    % Jacobian: at the fixed point (re,ri)
    Aee = -1/tauE + wEE*re*(1-re);
    Aei = -wEI*re*(1-re);
    Aie = wIE*ri*(1-ri);
    Aii = -1/tauI - wII*ri*(1-ri);
    
    A = [Aee Aei; Aie Aii];
    
    trA = trace(A);
    detA = det(A);
    
    d1 = 1/2*(trA + sqrt( trA^2 - 4*detA ));
    d2 = 1/2*(trA - sqrt( trA^2 - 4*detA ));
    
    Tr(k) = trA;
    De(k) = detA;
    out.Eig(k,j,1) = d1;
    out.Eig(k,j,2) = d2;

    r = [.6; .2];
    R = zeros(Tds,N);
    tt = 0;
    
    for t = 1:L    
      u = W*r + Io;
      K = feval(F,u);
      r = r + dt*(-r./tau + K);   
      if mod(t,ds)==0
        tt=tt+1;
        R(tt,:)=r;
      end 
    end
    
    RstatE = R(end-500:end,1);
    out.Re(k,j,1) = max(RstatE);
    out.Re(k,j,2) = min(RstatE);
    
    %         figure
    %         plot(R)
    %         set(gca,'ylim',[0 1])
    %         %title(sprintf('wEE = %g,  wEI = %g, v1 = %2.3f + i*%2.3f',wEE,wEI,real(d1),imag(d1)))
    %         title(sprintf('Ii = %2.2f, v1 = %2.3f + i*%2.3f',Ii,real(d1),imag(d1)))
    %         text(Tds/2,.3,sprintf('TrA = %2.3f, D = %2.3f',trA,trA^2-4*detA))
    
    if (Io(1) == -1.4 && Io(2) == -4) || (Io(1) == -1.4 && Io(2) == -4.8) || (Io(1) == -1.4 && Io(2) == -4.3)
      
      % null-clines (lines for which dr/dt=0)
      % 0 = -re/tauE + Fe( wEE*re - wEI*ri + Ie );
      % 0 = -ri/tauI + Fi( wIE*re - wII*ri + Ii );
      % FinvE( re/tauE ) = wEE*re - wEI*ri + Ie;
      % FinvI( ri/tauI ) = wIE*re - WII*ri + Ii;
      rE = 0:0.001:1;
      rI = 0:0.001:1;
      nullE(:,1) = ( wEE*rE + Io(1) - feval( FinvE, rE/tauE ) )/wEI; % = rI
      nullI(:,1) = ( wII*rI - Io(2) + feval( FinvI, rI/tauI ) )/wIE; % = rE
      
      trayectory(:,1) = R(1:10000,1) + 1i*R(1:10000,2);
      figure; set(gcf,'color','w'); plot(R(1:1000,1)); tp_editplots
      title(sprintf('Ii = %.3f',Io(2)))
      print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_bifurcation_timeseries_damped_j%d.pdf',j))
      
      figure_w;

      plot(rE,nullE(:,1),'k')
      hold on
      plot(nullI(:,1),rI,'color',[.5 .5 .5])
      plot(trayectory(:,1),'k:')
      % set(gca,'XTick',1:10:length(Iis ),'XTickLabels',num2cell(Iis(1:10:end)))
      axis(gca,'square')
      tp_editplots
      axis([0 0.8 0 0.8])
      
      print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_bifurcation_phasediag_j%d.pdf',j))


%     elseif Io(1) == -3
%       
%       rE = 0:0.001:1;
%       rI = 0:0.001:1;
%       nullE(:,2) = ( wEE*rE + Io(1) - feval( FinvE, rE/tauE ) )/wEI; % = rI
%       nullI(:,2) = ( wII*rI - Io(2) + feval( FinvI, rI/tauI ) )/wIE; % = rE
      
%       trayectory(:,2) = R(1:10000,1) + 1i*R(1:10000,2);
%       
%       figure; set(gcf,'color','w'); plot(R(1:1000,1)); tp_editplots
%       print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_bifurcation_timeseries_sustained.pdf'))
% %       
    end
  end
end

save(sprintf('~/pmod/proc/pmod_bifurcation_v%d.mat',v),'out')

error('!')

%%

load(sprintf('~/pmod/proc/pmod_bifurcation_v%d.mat',v))

close
k = 1;
figure; set(gcf,'color','white'); hold on
plot(out.Re(k,:,1),'k','linewidth',2)
plot(out.Re(k,:,2),'k','linewidth',2)
set(gca,'XTick',1:10:length(Iis),'XTickLabels',num2cell(Iis(1:10:end)))
xlabel('Inhibitory input'); ylabel('r_E')
% title(sprintf('Ies = %.3f',Ies(k)))
axis([1 121 0 0.4]); tp_editplots
line([find(Iis==-4) find(Iis==-4)],[0 0.4],'linestyle',':','color','k')
line([find(Iis==-4.8) find(Iis==-4.8)],[0 0.4],'linestyle',':','color','k')

print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_bifdiagram.pdf'));
%%

load(sprintf('~/pmod/proc/pmod_bifurcation_v%d.mat',1))

for i = 1 : size(out.Eig,1)
  for j = 1 : size(out.Eig,2)
    i
    % Stable node
    if all(isreal(out.Eig(i,j,:))) && all(out.Eig(i,j,:) < 0)
      out.re(i,j,1) = 1;
    % Stable focus: Noise driven
    elseif all(~isreal(out.Eig(i,j,:))) && all(real(out.Eig(i,j,:)) < 0)
      out.re(i,j,1) = 2;
    % Unstable focus: Oscillations?
    elseif all(~isreal(out.Eig(i,j,:))) && all(real(out.Eig(i,j,:)) > 0)
      out.re(i,j,1) = 3;
    % Unstable node
    elseif all(isreal(out.Eig(i,j,:))) && all(real(out.Eig(i,j,:)) > 0)
      out.re(i,j,1) = 4;
    % Saddle node
    elseif all(isreal(out.Eig(i,j,:))) && ~all(real(out.Eig(i,j,:)) > 0)
      out.re(i,j,1) = 5;
    else
      error('!')
    end
    
  end
end


figure; set(gcf,'color','w');hold on
subplot(1,1,1); %hold on
imagesc(out.re)
set(gca,'YTick',1:10:length(Iis ),'YTickLabels',num2cell(Iis(1:10:end)))
set(gca,'XTick',1:10:length(Ies),'XTickLabels',num2cell(Ies(1:10:end)))
axis(gca,'square')


%%


%%
figure
xSize = 15; ySize = 4.8;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[50 50 xSize*50 ySize*50],'Color','w')

axes('position',[.1 .2 .22 .69])
plot(Ies,out.Re(:,1),'k')
hold on
Ibif = Ies(find(abs(t(:,1)-out.Re(:,2))>0.0001,1,'first')); %bifurcations
plot(Ibif(1)*[1 1],[0 1],'k')
plot(Ibif(end)*[1 1],[0 1],'k:')
xlabel('Control param. \itI_{E}','fontname','times')
ylabel('\itr_E    ','fontname','times','rotation',0)
text(-.2,1.1,'A','units','normalized','fontsize',14)
tp_editplots
axes('position',[.42 .2 .22 .69])
plot(rE,nullE(:,1),'k')
hold on
plot(nullI(:,1),rI,'color',[.5 .5 .5])
plot(trayectory(:,1),'k:')
xlabel('\itr_E','fontname','times')
ylabel('\itr_I    ','rotation',0,'fontname','times')
set(gca,'xlim',[0 1],'ylim',[0 1],'xtick',0:.2:1,'ytick',0:.2:1)
title('\itI_{E}\rm = 0','fontname','times')
text(-.2,1.1,'B','units','normalized','fontsize',14)

tp_editplots
axes('position',[.74 .2 .22 .69])
plot(rE,nullE(:,2),'k')
hold on
plot(nullI(:,2),rI,'color',[.5 .5 .5])
plot(trayectory(:,2),'k:')
xlabel('\itr_E','fontname','times')
ylabel('\itr_I    ','rotation',0,'fontname','times')
set(gca,'xlim',[0 1],'ylim',[0 1],'xtick',0:.2:1,'ytick',0:.2:1)
title('\itI_{E}\rm = 2.5','fontname','times')
text(-.2,1.1,'C','units','normalized','fontsize',14)
tp_editplots
print(gcf,'-dpdf',sprintf('~/pmod/plots/pmod_bifurcation_v%d.pdf',1))%
%
%
%
%


