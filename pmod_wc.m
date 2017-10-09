% pmod_wc

clear all;
rng shuffle


% -------------------------------
% VERSION 1
% -------------------------------
v    = 2;
GAIN = 0:.02:1;
SIG  = 0.01;
Ies  = -3.4:0.1:0;
% -------------------------------
% VERSION 5
% -------------------------------
% v    = 5;
% GAIN = 0.02:0.02:0.5;
% SIG  = 0:0.005:0.02;
% Iis  = -4.4:0.2:-3.2;
% Ies  = -3:0.2:0;
% -------------------------------
% VERSION 6
% -------------------------------
% v    = 6;
% GAIN = 0.02:0.02:0.5;
% SIG  = 0.01;
% Iis  = -5:0.2:-3;
% -------------------------------
% GAIN = 0:.01:0.2;
% SIG  = 0.00:0.002:0.01;
p = 1.1;
% fixed params:

N = 4;

wII=4;
wIE=16;
wEI=10;
wEE=12;
g = 1;

tauE = 1;
tauI = 1;
tau = [tauE;tauI;tauE;tauI];

Iis = -4;

numIis = length(Iis);
numIes = length(Ies);
Io=zeros(N,1);

dt=0.05;
tmax = 1000;
tspan=0:dt:tmax;
L = length(tspan);

ds = 2;
Tds = length(0:ds*dt:tmax)-1;

% Connectivity:s
W11 = [wEE -wEI; wIE -wII];
W22 = W11;
W12 = [g 0; 0 0];
W21 = [g 0; 0 0];

W = [W11 W12; W21 W22];

nTrials = 10;

uu = -4:.1:4;
Transf = zeros(length(uu),2);

%%


for Q = 1:2  % conditions (no gain change, equivalent gain change, daE<daI, daE>daI)
  for k = 1:numIis
    for l = 1:numIes
      for ig = 1: length(GAIN)
        for ip = 1 : length(SIG)
          
          if ~exist(['~/pmod/proc/' sprintf('pmod_wc_q%d_ii%d_ie%d_g%d_p%d_v%d.txt',Q,k,l,ig,ip,v)])
            system(['touch ' '~/pmod/proc/' sprintf('pmod_wc_q%d_ii%d_ie%d_g%d_p%d_v%d.txt',Q,k,l,ig,ip,v)]);
          else
            continue
          end
          
          Rc     = 0;
          RcAmpl = 0;
          
          for trial = 1:nTrials
            fprintf('Starting computation trial %d ...\n',trial);

            Gain = GAIN(ig);
            %           p    = P(ip);
            sigma = SIG(ip);
            
            sdt = sqrt(dt)*sigma;
            
            Ie = Ies(l);
            Ii = Iis(k);
            Io(2) = Ii;
            Io(4) = Ii;
            Io(1) = Ie;
            Io(3) = Ie;
            
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
            elseif Q == 5
              gE = 1;
              gI = 1+p*Gain;
            elseif Q == 6
              gE = 1+p*Gain;
              gI = 1;
            end
            
            aE = 1/gE;
            aI = 1/gI;
            
            Fe = @(x) 1./(1 + exp(-x/aE) );
            Fi = @(x) 1./(1 + exp(-x/aI) );
            F = @(x) [feval(Fe,x(1));feval(Fi,x(2));feval(Fe,x(3));feval(Fi,x(4))];
            
            % store functions:
            Transf(:,1) = feval(Fe,uu);
            Transf(:,2) = feval(Fi,uu);
            
            if ~exist(sprintf('~/pmod/proc/pmod_wc_transfun_q%d_v%d.mat',Q,v))
              save(sprintf('~/pmod/proc/pmod_wc_transfun_q%d_v%d.mat',Q,v),'Transf');
            end
            
            r = rand(N,1);
            R = zeros(Tds,N);
            tt = 0;
            
            % Warm-up:
            for t = 1:2000
              u = W*r + Io;
              K = feval(F,u);
              r = r + dt*(-r./tau + K) + sdt*randn(N,1);
            end
            
            for t = 1:L
              u = W*r + Io;
              K = feval(F,u);
              r = r + dt*(-r./tau + K) + sdt*randn(N,1);
              if mod(t,ds)==0
                tt=tt+1;
                R(tt,:)=r;
              end
            end
            
            RstatE = R(:,[1 3]);
            Ampl(:,1)=abs(R(:,1) + 1i*R(:,2));
            Ampl(:,2)=abs(R(:,3) + 1i*R(:,4));
            
            rc = corr(RstatE);
            Rc = Rc+rc(2)/nTrials;
            
            rc = corr(RstatE);
            RcAmpl = RcAmpl+rc(2)/nTrials;
            
          end
          
          fprintf('Saving...\n');
          save(sprintf('~/pmod/proc/pmod_wc_q%d_ii%d_ie%d_g%d_p%d_v%d.mat',Q,k,l,ig,ip,v),'R','Rc','RcAmpl','-v7.3')
          
          try
            load(sprintf('~/pmod/proc/pmod_wc_q%d_ii%d_ie%d_g%d_p%d_v%d.mat',Q,k,l,ig,ip,v));
          catch me
            pause(randi(10))
            save(sprintf('~/pmod/proc/pmod_wc_q%d_ii%d_ie%d_g%d_p%d_v%d.mat',Q,k,l,ig,ip,v),'R','Rc','RcAmpl','-v7.3')
          end
          
        end
      end
    end
  end
end

error('!')

%%


clear Rc RcAmpl tmp_rc tmp_rcamp

% tmp_rc = zeros(26,2,length(GAIN),length(SIG));
% tmp_rcamp = zeros(26,2,length(GAIN),length(SIG));


for Q = 2:2  % conditions (no gain change, equivalent gain change, daE<daI, daE>daI)
  for k = 1:numIis
    for l = 1:numIes
      for ig = 1: length(GAIN)
        for ip = 1 : length(SIG)
           
%           try
            load(sprintf('~/pmod/proc/pmod_wc_q%d_ii%d_ie%d_g%d_p%d_v%d.mat',Q,k,l,ig,ip,v))
%           catch me
%             delete(sprintf('~/pmod/proc/pmod_wc_q%d_ii%d_g%d_p%d_v%d.txt',Q,k,ig,ip,v))
%             warning(sprintf('%s',me.message))
%           end

          tmp_rc(k,l,Q,ig,ip)    = Rc;
          tmp_rcamp(k,l,Q,ig,ip) = RcAmpl;
        
        end
      end
    end
  end
end

Rc_all = squeeze(tmp_rcamp(:,:,2,:,:));
% RcAmpl = tmp_rcamp;


%%
imagesc(Rc_all(:,:)',[-1 1])
set(gca,'xtick',[1:3:35],'xticklabel',Ies(1:3:35))
set(gca,'ytick',[1:3:51],'yticklabel',GAIN(1:3:51))


%%


Rc = squeeze(Rc_all(:,:,:,:));

addpath ~/Documents/MATLAB/cbrewer/cbrewer/
cmap = cbrewer('div', 'RdBu', 100,'pchip');% colormap(autumn)
cmap = cmap(end:-1:1,:);


diffs = [1 2; 1 3; 1 4; 3 4]
figure_white(10,10)

% xSize = 20; ySize = 20;
% xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
% set(gcf,'PaperUnits','centimeters')
% set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
% set(gcf,'Position',[50 50 xSize*50 ySize*50],'Color','w')

% idiff = 1: equal gain vs. no gain
% idiff = 2: E>I vs. no gain
% idiff = 3: I>E vs. no gain

for idiff = 1 : 1
  for k = 10 : 10%length(SIG)
    
%     subplot(5,4,k)
%     d = squeeze(Rc(k,diffs(idiff,2),:,2:end)-Rc(k,diffs(idiff,1),:,2:end));
        d = squeeze(Rc(:,diffs(idiff,2),:,k)-Rc(:,diffs(idiff,1),:,k));

    d=d.*(d>0.03 & squeeze(Rc(:,1,:,k))<0.8 & squeeze(Rc(:,2,:,k))>0.4 | d<-0.03 & squeeze(Rc(:,1,:,k))<0.8 & squeeze(Rc(:,2,:,k))>0.4);
    
%     imagesc(d(3:15,(3:15))',[-0.2 0.2]); colormap(cmap); axis square
        imagesc(d',[-0.2 0.2]); colormap(cmap); axis square

    title(['\itI_{i} : ' sprintf('%.2f',k)])
    
  end
%   print(gcf,'-djpeg100',sprintf('~/pmod/plots/pmod_wc_corrdiff_l%d_d%d_v%d.jpeg',l,idiff,v));
end

set(gca,'xtick',3:2:15,'xticklabel',Iis(3:2:15))
set(gca,'ytick',3:2:15,'yticklabel',GAIN(3:2:15))
set(gca,'tickdir','out')
xlabel('Bifurcation parameter (Inh. input)');
ylabel('Gain modulation');

print(gcf,'-depsc2',sprintf('~/pmod/plots/pmod_wc_corrdiff_nothresh_v%d.eps',v));

% Rc_all = Rc;



%%

p=14; g = 20 ;

Rcc = Rc_all(:,:,g,p);



figure

xSize = 15; ySize = 4.8;
xLeft = (21-xSize)/2; yTop = (30-ySize)/2;
set(gcf,'PaperUnits','centimeters')
set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
set(gcf,'Position',[50 50 xSize*50 ySize*50],'Color','w')

axes('position',[.1 .2 .22 .69])
plot([Iis(1) Iis(end)],[0 0],'k:')
hold on
plot(Iis,Rcc,'.-')
xlabel('Control param. \itI_{I}','fontname','times')
ylabel('corr. EE (c_{E-E })','fontname','times')
l = {'no gain','\deltaG_{E} = \deltaG_{I}',...
  '\deltaG_{E} > \deltaG_{I}','\deltaG_{E} < \deltaG_{I}'};
col=lines(4);
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
for i=1:4
  plot(xlim(1)+.46*(xlim(2)-xlim(1)),ylim(1)+(.38-(i-1)*.1)*(ylim(2)-ylim(1)),...
    'o','color',col(i,:),'markerfacecolor',col(i,:))
  text(.5,.38-(i-1)*.1,l{i},'units','normalized')
end
plot(-4.4*[1 1],ylim,'k-')
text(-.18,1.05,'A','units','normalized','fontsize',14)


Q = 3;
Iopt = find( Rcc(:,1)-Rcc(:,Q)>.05 & Rcc(:,1)<.7 & Rcc(:,Q)>0 );

if ~isempty(Iopt)
  plot(Iis(Iopt)'*[1 1],ylim,'k:')
end


load(sprintf('~/pmod/proc/pmod_wc_transfun_q%d_v%d.mat',1,v))

axes('position',[.39 .6 .22 .3])
plot(uu,Transf(:,1),'b')
hold on

load(sprintf('~/pmod/proc/pmod_wc_transfun_q%d_v%d.mat',Q,v))

plot(uu,Transf(:,1),'k')
ylabel('F_{E}','fontname','times')
text(-.18,1.1,'B','units','normalized','fontsize',14)

load(sprintf('~/pmod/proc/pmod_wc_transfun_q%d_v%d.mat',1,v))

axes('position',[.39 .2 .22 .3])
plot(uu,Transf(:,2),'b')
hold on

load(sprintf('~/pmod/proc/pmod_wc_transfun_q%d_v%d.mat',Q,v))

plot(uu,Transf(:,2),'k')
xlabel('Input','fontname','times')
ylabel('F_{I}','fontname','times')


aE = 1;
aI = 1;

Fe = @(x) 1./(1 + exp(-x/aE) );
Fi = @(x) 1./(1 + exp(-x/aI) );
F = @(x) [feval(Fe,x(1));feval(Fi,x(2));feval(Fe,x(3));feval(Fi,x(4))];

Ii = Iis(Iopt(1));
Io(2) = Ii;
Io(4) = Ii;

r = rand(N,1);
R = zeros(Tds,N);
tt = 0;
for t = 1:floor(L/2)
  u = W*r + Io;
  K = feval(F,u);
  r = r + dt*(-r./tau + K) + sdt*randn(N,1);
  if mod(t,ds)==0
    tt=tt+1;
    R(tt,:)=r;
  end
end

RstatE = R(tt-1999:tt,[1 3]);

axes('position',[.68 .6 .3 .3])

time = ds*dt*(0:size(RstatE,1)-1);
plot(time,RstatE)
text(-.18,1.1,'C','units','normalized','fontsize',14)

aE = 1/(1+p*Gain);
aI = 1/(1+Gain);

Fe = @(x) 1./(1 + exp(-x/aE) );
Fi = @(x) 1./(1 + exp(-x/aI) );
F = @(x) [feval(Fe,x(1));feval(Fi,x(2));feval(Fe,x(3));feval(Fi,x(4))];

Ii = Iis(Iopt(1));
Io(2) = Ii;
Io(4) = Ii;

r = rand(N,1);
R = zeros(Tds,N);
tt = 0;
for t = 1:floor(L/2)
  u = W*r + Io;
  K = feval(F,u);
  r = r + dt*(-r./tau + K) + sdt*randn(N,1);
  if mod(t,ds)==0
    tt=tt+1;
    R(tt,:)=r;
  end
end

RstatE = R(tt-1999:tt,[1 3]);

axes('position',[.68 .2 .3 .3])

time = ds*dt*(0:size(RstatE,1)-1);
plot(time,RstatE)

text(.4,.7,sprintf('I_{I} = %2.2f',Ii),'units','normalized')

% export_fig Fig_TwoWC_Correlations.tiff -r300