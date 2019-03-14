
%-------------------------------------------------------------------------
% VERSION 1: 20-10-2018
% %-------------------------------------------------------------------------
v           = 1;
Ies         = -4:0.025:-1;
Iis         = -5:0.025:-2;
Gg          = 0:0.1:3;
Gains       = 0; 
nTrials     = 1;
tmax        = 3750;  % in units of tauE
EC          = 0;
%-------------------------------------------------------------------------
% VERSION 2: 20-10-2018: DETERMINE GLOBAL COUPLING PARAMETER
% %------------------------------------------------------------------------
% v           = 2;
% Ies         = -4:0.025:-1;
% Iis         = -5:0.025:-2;
% Gg          = 1.7;
% Gains       = [0 0.025:0.025:0.4 -0.025:-0.025:-0.4]; 
% nTrials     = 1;
% tmax        = 3750;  % in units of tauE
% EC          = 0;
%--------------------------------------------------------------------------

outdir = '~/pmod/proc/';

vv = v;

for iies = 1 : length(Ies)
  for iiis = 1 : length(Iis)
    for iG = 1:length(Gg)
      for igain = 1:length(Gains)
        
        load(sprintf([outdir 'pmod_wc_wholebrain_detosc_Ie%d_Ii%d_G%d_gain%d_v%d.mat'],iies,iiis,iG,igain,vv))
        osc1(iies,iiis,iG,igain) = mean(squeeze(mean(squeeze(out.osc1),1)));
        
      end
    end
  end
end
  
save(sprintf([outdir 'pmod_wc_wholebrain_detosc_all_v%d.mat'],vv),'osc1')
  
% DELETE OLD FILES
for iies = 1 : length(Ies)
  iies
  for iiis = 1 : length(Iis)
    for iG = 1:length(Gg)
      for igain = 1:length(Gains)
        
        warning('Deleting...')
        delete(sprintf([outdir 'pmod_wc_wholebrain_detosc_Ie%d_Ii%d_G%d_gain%d_v%d.mat'],iies,iiis,iG,igain,vv))
        delete(sprintf([outdir 'pmod_wc_wholebrain_detosc_Ie%d_Ii%d_G%d_gain%d_v%d_processing.txt'],iies,iiis,iG,igain,vv))
      end
    end
  end
end
