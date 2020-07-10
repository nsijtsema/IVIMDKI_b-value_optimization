% This scripts finds optimal DWI b-values for given expected parameter
% value-ranges, noise leveland provided model based on the Cramér-Rao Lower
% bound. A correction for echo time (TE) is applied to correct for minimum
% achievable TE at varying maximum b-values.
%
% Copyright (C)2020 Nienke Sijtsema
% Contact: n.sijtsema@erasmusmc.nl
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

%% Specify set of parameters:
nin = 36*16; % number of test parameters

ADC_range=[0.25 3.41]*1e-3;
perf_frac_range=[0.09 0.42];
Dstar_range= [6.29 237.39]*1e-3;  %decrease to 237.39 or less (10% left)    20 of signal left at b=10 for D*=0.16                                                                                                               
S0_range=[1500 2500]*(1/0.350813688); %Increase to 1500-2500, Correction factor for T2 decay.
K_range = [0.1 2.81];  %0.5 3 %2.81

sigma =100; %NL20 for paper!!! 1/15;  % specify noise level (for CRLB analysis)

prims = primes(11); % all prime number uptil and including input
ADC = HaltonSequence( nin, prims(1), .5) * (ADC_range(2)-ADC_range(1))+ADC_range(1); %HaltonSequence outputs pseudorandom numbers between 0 and 1. ADC is a scaled in the range given above.
perf_frac = HaltonSequence( nin, prims(2), .5) * (perf_frac_range(2)-perf_frac_range(1))+perf_frac_range(1);
Dstar = HaltonSequence( nin, prims(3), .5) * (Dstar_range(2)-Dstar_range(1))+Dstar_range(1);
K = HaltonSequence( nin, prims(4), .5) * (K_range(2)-K_range(1))+K_range(1);
S0= HaltonSequence( nin, prims(5), .5) * (S0_range(2)-S0_range(1))+S0_range(1);

xtst = [S0;perf_frac;Dstar;ADC;K]; 
predfun = @pred_NGIVIM_fNDS; %Choose prediction model -> is input for CRLB
varnames = { 'S0','f','Dstar','ADC','K'};

%% specify settings range
minb = 0;  % minimum b-value
maxb = 1500; % maximum b-value

minsett = minb;
maxsett = maxb;


initsett = [0 10 20 40 60 80 100 120 150 170 220 230 410 420 610 620 790 800 940 950 960 1200 1220 1230 1430 1440 1450 1460 1470 1480]';

% Exclude xtst samples where ADC & Dstar are too close together to separate
% (Dstar<6*ADC) and exclude nonphysical result where the curve goes upward
% (derivative>0)

delete_index= Dstar<(6*ADC);

dpred_db_maxb=S0.*((1-perf_frac).*exp(-maxb .* ADC + (1/6) .* maxb.^2 .* ADC.^2 .* K).*((1/3)*maxb.*ADC.^2.*K-ADC)-perf_frac.*Dstar.*exp(-maxb .* Dstar)); %Checked with wolfram alpha
delete_index = delete_index | dpred_db_maxb>0 ;

xtst(:,delete_index)=[];

whos xtst

%% Create optimization function:

% Adjust of each image the b-value. 

% compute CRLB of initial settings, only used to create selector
[CRLB,I,J] = CramerRaoLowerBound_MRI( xtst(:,1), @(x) predfun(x, initsett), sigma);  %inputs parameters, function and noise level. x corresponds with xtxt(:,1) (?)

varidx = [2 3 4 5];% we want to measure the CRLB of f, Dstar, D and K

selector = zeros(numel(varidx ),numel(I));
for k=1:numel(varidx )
    selector(k,I==J & I==varidx(k) ) = 1; 
end;

%all function handles! (needed for cost function: incl CRLB, acqtime and
%constraint for negative b-values.
TE= @(sett) 0.015*max(sett)+61; %TE in ms.
T2DecayAtTE= @(sett) exp(-TE(sett)./80);
ApplyT2DecayToTht= @(sett) [xtst(1,:)*T2DecayAtTE(sett); xtst(2:end,:)];  %S0 f D* D K
absCRLB = @(sett) (selector * CramerRaoLowerBound_MRI( ApplyT2DecayToTht(sett), @(x) predfun(x, sett), sigma )); 
relCRLB = @(sett) (selector * CramerRaoLowerBound_MRI( ApplyT2DecayToTht(sett), @(x) predfun(x, sett), sigma ))./xtst(varidx,:).^2; %CRLB relative to original values.
acqtime = @(sett) (TE(sett)+120)*numel(sett); %Acquisition time = (TE + "Readout time")*amount of b values.  RO=TR/num_slices-TE
constrain = @(sett) -sum( sett(sett<minb)-minb)+sum( sett(sett>maxb)-maxb) ; % add (severe) cost for negative b values and larger than 1500. 
relCRLBweighted=@(sett) relCRLB(sett).*[1;1;1;1];
costfun = @(sett) sqrt( mean(mean( relCRLBweighted(sett).^2 ,1),2) ) * acqtime(sett) + 1000*constrain( sett );


%% Perform local search from initial settings:

tst_sett = fminsearch( costfun, initsett(:)); 

%% plot contribution to relative CRLB for each parameter for all test points:
figure
set(gcf,'Name','testset');
plot( sqrt( relCRLBweighted( tst_sett)') ) %was tst_sett
legend( varnames{varidx} ) 
title('relative standard deviation')


figure
set(gcf, 'Name','testset');
plotselidx = 40;%sort(ceil(rand(1,5)*size(xtst,2))); %linspace(1,size(xtst,2),4));
plotselx = xtst(:, plotselidx) ; % select a few quasi random elements from xtst. 
plot_b = linspace(0,sqrt(2000),100).^2'; % use square to get more samples at low b-values. 
lh1 = plot( tst_sett , predfun( plotselx, tst_sett ) ,'*' ); 
hold on
lh2 = plot( plot_b , predfun( plotselx, plot_b) ,'-' ); 
hold off
legstr  = cell(1, numel(plotselidx ));
for k= 1 : numel(plotselidx )
    legstr{k} =  sprintf('test voxel %d',plotselidx (k) );
end
legend(lh1,legstr{:})
title(sprintf('Predicted signal intensity for %d selected voxels.',numel(plotselidx)))
xlabel('b-value');
ylabel('S (a.u)');

