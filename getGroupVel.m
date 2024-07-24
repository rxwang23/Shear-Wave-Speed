function getGroupVel(arfidata,axial,lat,t, iteration, C1_x)

indexFocus = knnsearch(lat,0);
% indexFocus = knnsearch(lat,24.5);

aimg.swPushLoc = indexFocus;

aimg.focalDepthmm = knnsearch(axial,0);
aimg.PRFKHz = 1/mean(diff(t)); 
aimg.numPreFrames = 1;

filePath = pwd;

saveFinal = 1;
blat=lat;
% factor.dx=0.2e-3;
factor.dx=0.5e-3;
factor.spaDs=1;
factor.NT=4096;
factor.axAvg=2e-3; % [m] axial average ROI in 
factor.tAvg=7; % [sample] axial average ROI in t
factor.tWind1=0.3; % tuckey window value
factor.tWind2=0.3; % tuckey window value
factor.tInterFac=2;
factor.NT=round(4096*factor.tInterFac/2);
factor.Dir_cutoffsVel=[30 0.2]; % in m/s
factor.maxFreq=2e3;
factor.Dir_order=2; % bandpass order
factor.Dir_power=2; % bandpass power
factor.axW=4.0e-3;factor.laW=factor.axW;
factor.windSize=factor.axW;
factor.windSizeCC=round(factor.windSize/factor.dx);
factor.stepSize=round(factor.axW/factor.dx);
for uu=1
    pushL=aimg.swPushLoc(uu);
     % if pushL<64
         factor.Dir_angles=0;
         direc='lr';
    % else
    %      factor.Dir_angles=pi;
    %      direc='rl';
    %  end

    factor.pushLoc=blat(pushL)*1e-3;
    factor.latCheck=round(pushL * 0.2); % this needs to be max 300
    disp("rmRevIsoPixel start");
    [data,axial,lat,t]=rmRevIsoPixel(arfidata,axial,lat,t,factor,aimg);
    disp("rmRevIsoPixel done");

    disp(2*[50 200]*median(diff(t)));
        disp(median(diff(t)));

    [b,a]=butter(2,2*[50 200]*median(diff(t)),'bandpass'); 
    
    dataF=permute(filtfilt(b,a,permute(data,[3 1 2])),[2 3 1]);
    dataF=movmean(dataF,round(factor.axAvg/factor.dx)); % axial direction
    dataF=permute(movmean(permute(dataF,[3 1 2]),factor.tAvg),[2 3 1]); % 5 point along time
    % dataF=diff(dataF,1,3); % particle velocity
    % t=t(1:end-1);
    tWin=tukeywin(size(dataF,3),factor.tWind1);
    tWinMat=permute(repmat(tWin,[1 size(dataF(:,:,1))]),[2 3 1]);
    dataw=dataF.*tWinMat;
    fixedHiFreq=factor.maxFreq;
    factor.Dir_cutoffs=fixedHiFreq*factor.dx*0.5./factor.Dir_cutoffsVel;
    disp("df2s start");

    dataw = df2d_Song_V4(dataw, factor.dx, factor.dx, factor.Dir_cutoffs,factor.Dir_order,factor.Dir_power,factor.Dir_angles);
    disp("df2s end");

    prf=1/mean(diff(t));          
    disp("sws start");
    [cs,cc,csx,csy] = CUSEShearWaveSpeed2DCCAHParFor(dataw,prf,factor.dx,factor.dx,direc,factor.stepSize,factor.windSizeCC,factor.tInterFac);
    disp("sws end");

    if saveFinal
        saveName = strcat(num2str(iteration), "_C1_", "C2_",  "R_", "eta");
        if ~exist(fullfile(filePath,'3D_200_homogenous'))
             mkdir(fullfile(filePath,'3D_200_homogenous'))
        end
    delete(fullfile(filePath,'3D_200_homogenous',saveName));
    save(fullfile(filePath,'3D_200_homogenous',saveName),'dataw','lat','axial','t','cs','csx','csy','cc','factor')
    end
end

filename = fullfile(filePath,'3D_200_homogenous',saveName);

figure('Name','SMS figure');
imagesc(lat,axial,csx);
hold on;
viscircles([C1_x,0], 5e-3,'Color','r');
colormap('jet'); colorbar;
xlabel('Laterial [m]');
ylabel('Axial [m]');
title('cs-x');
clim([0 3]);

xlim([0 0.04]);
ylim([-0.02 0.02]);

set(gca,'fontsize',20); set(gca,'fontweight','bold');
saveas(gcf,strcat(filename,"-cs-x.png"));

return