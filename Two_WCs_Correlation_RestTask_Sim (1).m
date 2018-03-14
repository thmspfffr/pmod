
% Stochastic simulation of 2 WC nodes subject to drug (or not) in Rest and
% Task

%-------------------------------------------------------------------------


clear all;

Drug = 2;
nTrials = 100;

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

dt=0.01;
tmax = 80000/400; % in units of tauE
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

% Control params.
%--------------------
Ie = -2.85; 
Ii = -3.50; %-1;

% input due to task:
dIe = 1.0;
dIi = 1.3;

% Drug effects:
if Drug==1
% inputs    
dIe_drug = 0;%0.08;
dIi_drug = -.08;%0;    
% gain modulation:
Gain_E = .20;
Gain_I = .20;
else
% inputs
dIe_drug = -0.20;
dIi_drug = -0.18;
% gain modulation:
Gain_E = .10; %0.07;
Gain_I = .10; %0.07;
end
    

Rate = zeros(nTrials,4);
%RateSD = zeros(4,1);
Cee = zeros(nTrials,4);
%CeeSD = zeros(4,1);

display('REST ...')     
%--------------------------------------------------------------------------

% transfer function:
gE = 1;
gI = 1;
aE = 1/gE;
aI = 1/gI;
Fe = @(x) 1./(1 + exp(-x/aE) );
Fi = @(x) 1./(1 + exp(-x/aI) );
F = @(x) [feval(Fe,x(1));feval(Fi,x(2));feval(Fe,x(3));feval(Fi,x(4))];

% Working point:
Io(1) = Ie;
Io(3) = Ie;        
Io(2) = Ii;
Io(4) = Ii;

c = zeros(1,nTrials);
for tr=1:nTrials
R = zeros(Tds,N);    
r = 0.001*rand(N,1);
tt = 0;
% transient:
for t = 1:10000
    u = W*r + Io;            
    K = feval(F,u);            
    r = r + dt*(-r + K)./tau + sqrt(dt)*sigma*randn(N,1);
end     
% simulation
 for t = 1:L            
     u = W*r + Io;            
     K = feval(F,u);            
     r = r + dt*(-r+ K)./tau + sqrt(dt)*sigma*randn(N,1);            
     if mod(t,ds)==0
        tt=tt+1;
        R(tt,:)  = r;
     end       
 end

%correlation: 
rE = R(:,[1 3]);
rc = corr(rE);
Cee(tr,1) = rc(2);
%CeeSD(1) = high(2)-low(2);

rEo = mean(rE(:,1));
%rEsd = std(rE(:,1));
Rate(tr,1) = rEo-rEo;
%RateSD(1) = rEsd;

end

display('REST + drug ...')     
%--------------------------------------------------------------------------

% transfer function:
gE = 1 + Gain_E;
gI = 1 + Gain_I;
aE = 1/gE;
aI = 1/gI;
Fe = @(x) 1./(1 + exp(-x/aE) );
Fi = @(x) 1./(1 + exp(-x/aI) );
F = @(x) [feval(Fe,x(1));feval(Fi,x(2));feval(Fe,x(3));feval(Fi,x(4))];

% Working point:
Io(1) = Ie + dIe_drug;
Io(3) = Ie + dIe_drug;        
Io(2) = Ii + dIi_drug;
Io(4) = Ii + dIi_drug;

for tr=1:nTrials
r = 0.001*rand(N,1);
R = zeros(Tds,N);
tt = 0;
% transient:
for t = 1:10000
    u = W*r + Io;            
    K = feval(F,u);            
    r = r + dt*(-r + K)./tau + sqrt(dt)*sigma*randn(N,1);
end     
% simulation
 for t = 1:L            
     u = W*r + Io;            
     K = feval(F,u);            
     r = r + dt*(-r+ K)./tau + sqrt(dt)*sigma*randn(N,1);            
     if mod(t,ds)==0
        tt=tt+1;
        R(tt,:)  = r;
     end       
 end

%correlation: 
rE = R(:,[1 3]);
rc = corr(rE);
Cee(tr,2) = rc(2);
%CeeSD(1) = high(2)-low(2);

rEm = mean(rE(:,1));
%rEsd = std(rE(:,1));
Rate(tr,2) = rEm-rEo;
%RateSD(1) = rEsd;
end

display('TASK ...')     
%--------------------------------------------------------------------------

% transfer function:
gE = 1;
gI = 1;
aE = 1/gE;
aI = 1/gI;
Fe = @(x) 1./(1 + exp(-x/aE) );
Fi = @(x) 1./(1 + exp(-x/aI) );
F = @(x) [feval(Fe,x(1));feval(Fi,x(2));feval(Fe,x(3));feval(Fi,x(4))];

% Working point:
Io(1) = Ie + dIe;
Io(3) = Ie + dIe;        
Io(2) = Ii + dIi;
Io(4) = Ii + dIi;

for tr=1:nTrials
r = 0.001*rand(N,1);
R = zeros(Tds,N);
tt = 0;
% transient:
for t = 1:10000
    u = W*r + Io;            
    K = feval(F,u);            
    r = r + dt*(-r + K)./tau + sqrt(dt)*sigma*randn(N,1);
end     
% simulation
 for t = 1:L            
     u = W*r + Io;            
     K = feval(F,u);            
     r = r + dt*(-r+ K)./tau + sqrt(dt)*sigma*randn(N,1);            
     if mod(t,ds)==0
        tt=tt+1;
        R(tt,:)  = r;
     end       
 end

%correlation: 
rE = R(:,[1 3]);
rc = corr(rE);
Cee(tr,3) = rc(2);
%CeeSD(1) = high(2)-low(2);

rEm = mean(rE(:,1));
%rEsd = std(rE(:,1));
Rate(tr,3) = rEm-rEo;
%RateSD(1) = rEsd;
end

display('TASK + drug ...')     
%--------------------------------------------------------------------------

% transfer function:
gE = 1 + Gain_E;
gI = 1 + Gain_I;
aE = 1/gE;
aI = 1/gI;
Fe = @(x) 1./(1 + exp(-x/aE) );
Fi = @(x) 1./(1 + exp(-x/aI) );
F = @(x) [feval(Fe,x(1));feval(Fi,x(2));feval(Fe,x(3));feval(Fi,x(4))];

% Working point:
Io(1) = Ie + dIe_drug + dIe;
Io(3) = Ie + dIe_drug + dIe;        
Io(2) = Ii + dIi_drug + dIi;
Io(4) = Ii + dIi_drug + dIi;

for tr=1:nTrials
r = 0.001*rand(N,1);
R = zeros(Tds,N);
tt = 0;
% transient:
for t = 1:10000
    u = W*r + Io;            
    K = feval(F,u);            
    r = r + dt*(-r + K)./tau + sqrt(dt)*sigma*randn(N,1);
end     
% simulation
 for t = 1:L            
     u = W*r + Io;            
     K = feval(F,u);            
     r = r + dt*(-r+ K)./tau + sqrt(dt)*sigma*randn(N,1);            
     if mod(t,ds)==0
        tt=tt+1;
        R(tt,:)  = r;
     end       
 end

%correlation: 
rE = R(:,[1 3]);
rc = corr(rE);
Cee(tr,4) = rc(2);
%CeeSD(1) = high(2)-low(2);

rEm = mean(rE(:,1));
%rEsd = std(rE(:,1));
Rate(tr,4) = rEm-rEo;
%RateSD(1) = rEsd;
end

%--------------------------------------------------------------------------    

load  WCs_Correlations_ParamSpace2 CeE p Gain W g tau Iis Ies

l = {'no modulation','\deltaG_{E} = \deltaG_{I}',...
    '\deltaG_{E} > \deltaG_{I}','\deltaG_{E} < \deltaG_{I}'};

figure
xSize = 12.5; ySize = 9;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[50 50 xSize*50 ySize*50],'Color','w')

axes('position',[.1 .63 .28 .28])
pos=get(gca,'position');
imagescnan(Iis,Ies,CeE(:,:,1),[0 .3])
axis xy
set(gca,'fontsize',8)
xlabel('\itI_{I}','fontname','times','fontsize',11)
ylabel('\itI_{E}  ','fontname','times','rotation',0,'fontsize',11)
h=colorbar;
set(h,'position',[pos(1)+pos(3)+.01 pos(2) .01 pos(4)],'fontsize',8)
box on
set(gca,'linewidth',1)
text(.23,1.13,'no modulation','units','normalized','fontsize',11)
text(-8,-1,'oscillations','fontsize',8)
text(1.25,.5,'C_{EE}','units','normalized','fontsize',9)
text(-.2,1.15,'A','units','normalized','fontsize',14)
set(gca,'xlim',[-12 0],'ylim',[-7 7],'xtick',-12:2:0,'ytick',-6:2:6)
hold on
plot(Ii,Ie,'rs')
plot(Ii+dIi,Ie+dIe,'ws')
plot(Ii+dIi_drug,Ie+dIe_drug,'rs','markerfacecolor','w')
plot(Ii+dIi+dIi_drug,Ie+dIe+dIe_drug,'ws','markerfacecolor','w')

axes('position',[.58 .63 .28 .28])
pos=get(gca,'position');
imagescnan(Iis,Ies,CeE(:,:,2)-CeE(:,:,1),[-.1 .1])
axis xy
set(gca,'fontsize',8)
xlabel('\itI_{I}','fontname','times','fontsize',11)
ylabel('\itI_{E}  ','fontname','times','rotation',0,'fontsize',11)
h=colorbar;
set(h,'position',[pos(1)+pos(3)+.01 pos(2) .01 pos(4)],'fontsize',8)
box on
set(gca,'linewidth',1)
text(-.02,1.13,['gain modulation ' '(' l{2} ')'],'units','normalized','fontsize',11)
text(-8,-1,'oscillations','fontsize',8)
text(1.25,.5,'\DeltaC_{EE}','units','normalized','fontsize',9)
text(-.2,1.15,'B','units','normalized','fontsize',14)
set(gca,'xlim',[-12 0],'ylim',[-7 7],'xtick',-12:2:0,'ytick',-6:2:6)
hold on
plot(Ii,Ie,'rs')
plot(Ii+dIi,Ie+dIe,'ws')
plot(Ii+dIi_drug,Ie+dIe_drug,'rs','markerfacecolor','w')
plot(Ii+dIi+dIi_drug,Ie+dIe+dIe_drug,'ws','markerfacecolor','w')

axes('position',[.1 .29 .35 .17])
rate = mean(Rate);
rateSD = std(Rate)/sqrt(nTrials);
errorbar(1:4,rate,rateSD,'k.')
hold on
bar(1,rate(1),'facecolor','b','barwidth',0.9)      
bar(2,rate(2),'facecolor',[.6 .6 .6],'barwidth',0.9)   
bar(3,rate(3),'facecolor',[.5 .5 1],'barwidth',0.9)        
bar(4,rate(4),'facecolor',[.9 .9 .9],'barwidth',0.9)       
text(-.27,.5,'\Deltar_{E}','units','normalized','fontsize',9)
box off
text(-.2,1.15,'C','units','normalized','fontsize',14)


axes('position',[.58 .29 .32 .17])
Corr = mean(Cee);
CorrSD = std(Cee)/sqrt(nTrials);
errorbar(1:4,Corr,CorrSD,'k.')
hold on
bar(1,Corr(1),'facecolor','b','barwidth',0.9)      
bar(2,Corr(2),'facecolor',[.6 .6 .6],'barwidth',0.9)   
bar(3,Corr(3),'facecolor',[.5 .5 1],'barwidth',0.9)        
bar(4,Corr(4),'facecolor',[.9 .9 .9],'barwidth',0.9)       

p12 = ranksum(Cee(:,1),Cee(:,2));
if p12<0.05, plot(2,Corr(2)+CorrSD(2)+.03,'k*'), end
p34 = ranksum(Cee(:,3),Cee(:,4));
if p34<0.05, plot(4,Corr(4)+CorrSD(4)+.03,'k*'), end

set(gca,'fontsize',8,'xtick',1:3,'ylim',[0 .15],'xtick',[])    
text(1,-.03,'Rest','fontsize',8,'horizontalalignment','center','color','b')
text(2,-.03,'Rest+drug','fontsize',8,'horizontalalignment','center')
text(3,-.03,'Task','fontsize',8,'horizontalalignment','center','color',[.5 .5 1])
text(4,-.03,'Task+drug','fontsize',8,'horizontalalignment','center','color',[.5 .5 .5])

text(-.27,.5,'C_{EE}','units','normalized','fontsize',9)
box off
text(-.2,1.15,'D','units','normalized','fontsize',14)
