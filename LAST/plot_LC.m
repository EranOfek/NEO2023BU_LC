%% plot_LC for 2023BU


%% load all MatchedSources files and calibrate photometry
F  = dir('MS*.mat');
Nf = numel(F);
AllLC = zeros(0,7);
for If=1:1:Nf
    MS = io.files.load2(F(If).name);

    [MSz, LC] = lcUtil.zp_external(MS);

    AllLC = [AllLC; LC];
end

Ephem = io.files.load2('Ephem.mat');

% AllLC columns:
ColAllLC = {'JD', 'FLUX_PSF', 'MAGERR_PS', 'Chi2dof', 'FLAGS', 'FLAG_POS', 'MAG_PSF', 'ErrZP'};

save -v7.3 AllLC.mat AllLC ColAllLC

%%
load AllLC.mat

RAD = 180./pi;
AU = 149597870.700;  % km
C  = 299792.458;
LT = AU/C;
SEC_DAY = 86400;


JD0 = 2459970;

% sort LC
AllLC = sortrows(AllLC,1);
FlagPos = AllLC(:,6)==1;
AllLC = AllLC(FlagPos,:);


JD  = AllLC(:,1);
AllLC(:,1) = AllLC(:,1) - JD0;

% clean the data

FlagGood = (AllLC(:,1)<1 & AllLC(:,2)<3100) | (AllLC(:,1)>1 & AllLC(:,2)<45000);
AllLC    = AllLC(FlagGood,:);
JD       = JD(FlagGood);


% estimate angular speed of asteroid
Ne = sizeCatalog(Ephem);
Dist = zeros(Ne,1);
for Ie=1:1:Ne-1
    Dist(Ie) = celestial.coo.sphere_dist_fast(Ephem.Catalog(Ie,2), Ephem.Catalog(Ie,3), Ephem.Catalog(Ie+1,2), Ephem.Catalog(Ie+1,3));
end
AngSpeed = Dist.*RAD.*3600./60;  % [arcsec/sec]

% ang dist to first point (ref)
Dist1 = celestial.coo.sphere_dist_fast(Ephem.Catalog(1,2), Ephem.Catalog(1,3),Ephem.Catalog(:,2), Ephem.Catalog(:,3));

% correct time for light time correction
Delta = interp1(Ephem.Catalog(:,1), Ephem.Catalog(:,4), JD);
r     = interp1(Ephem.Catalog(:,1), Ephem.Catalog(:,5), JD);
STO   = interp1(Ephem.Catalog(:,1), Ephem.Catalog(:,6), JD);
AngSpeed = interp1(Ephem.Catalog(:,1), AngSpeed, JD);
AngDist1 = interp1(Ephem.Catalog(:,1), Dist1, JD);

LightTimeCorr = Delta.*LT;

Period = 38.*2;
JD_AstRefFrame = JD - LightTimeCorr./SEC_DAY; % + Period.*AngDist1./(2.*pi)./SEC_DAY;

% coorect magnitude to distance
% to avoid excessive changes
Flux = AllLC(:,2);
FluxDistCorr = Flux.*Delta.^2; %./AngSpeed;

plot(JD_AstRefFrame, FluxDistCorr,'k-')

% average LC in order to subntrcat mean
B = timeseries.binning([JD_AstRefFrame, FluxDistCorr],20./1440,[NaN NaN],{'MidBin', @median, @std, @numel});
MeanInterp = interp1(B(:,1),B(:,2),JD_AstRefFrame);

FluxBiasRemoved = FluxDistCorr./MeanInterp;


Freq = (0:1e-6:0.1).';
[P,F] = timeseries.period([(JD_AstRefFrame-JD0).*SEC_DAY, FluxDistCorr],Freq);
Flag1 = (JD_AstRefFrame-JD0)<1;
Flag2 = (JD_AstRefFrame-JD0)>1;

[P1,F1] = timeseries.period([(JD_AstRefFrame(Flag1)-JD0).*SEC_DAY, FluxDistCorr(Flag1)],Freq);
[P2,F2] = timeseries.period([(JD_AstRefFrame(Flag2)-JD0).*SEC_DAY, FluxDistCorr(Flag2)],Freq);

F2.Per(end)-F1.Per(end)

%%

MagLC=LC(:,[1 7 3]);
F=MagLC(:,2)<19;
MagLC=MagLC(F,:);
MagLC = sortrows(MagLC,1);

MagLC(:,1) = (MagLC(:,1)-2459971.0).*1440;

plot.errorxy(MagLC,'FaceColor','k','EdgeColor','k')
hold on
plot(MagLC(:,1),MagLC(:,2),'k-','Color',[0.8 0.8 0.8])    

plot.invy;
H = xlabel('Time [min]');
H.FontSize = 18;
H.Interpreter = 'latex';
H = ylabel('$B_{\rm p}$ [mag]');
H.FontSize = 18;
H.Interpreter = 'latex';



