function [data,axial,lat,t]=rmRevIsoPixel(data,axial,lat,t,factor,aimg)
    PRF_interval=median(diff(t));
    diffT=round(diff(t)*100)/100;
    ind_1=[true ~(diffT~=(round((100*PRF_interval)))/100)]';
    ind_2=circshift(ind_1,1);
    ind_3=circshift(ind_2,1);
    ind=logical(ind_1.*ind_2.*ind_3);
    t_1=t(ind);
    data_1=data(:,:,ind);
    meanNoise=median(data_1(:,:,1:aimg.numPreFrames),3);
    meanNoiseRep=repmat(meanNoise,[1 1 size(data_1,3)]);
    data_1=data_1-meanNoiseRep;
    t_1=t_1-t_1(aimg.numPreFrames);
    t_interp_1=t_1(1):1/aimg.PRFKHz:t_1(end);
    FD=aimg.focalDepthmm-2; % reason?
    FD_ind=knnsearch(axial,FD);
    
    
    latId_1=factor.latCheck;
    
    dispThresh=10;
    % latId_1 = min(latId_1, size(data_1, 2));

    disp_temp = (squeeze(data_1(FD_ind,latId_1,:)));
    
    ind_1=[true;~(abs(diff(disp_temp))>dispThresh)]; % find interframe displacement > 1 micron
    if all(ind_1)==1
        try data_1=interp1(t_1,permute(data_1,[3 1 2]),t_interp_1,'spline');
        catch ME        
            for ii=1:size(data_1,1)
                for jj=1:size(data_1,2)
                    if all(isnan(squeeze(data_1(ii,jj,:))))==1
                     data_1(ii,jj,[1 size(data_1,3)])=[0 1];
                    end
                end
            end
            data_1=interp1(t_1,permute(data_1,[3 1 2]),t_interp_1,'spline');
        end
    else
        t_1=t_1(ind_1);
        data_1=data_1(:,:,ind_1);
        try data_1=interp1(t_1,permute(data_1,[3 1 2]),t_interp_1,'spline');
        catch ME
            for ii=1:size(data_1,1)
                for jj=1:size(data_1,2)
                    if all(isnan(squeeze(data_1(ii,jj,:))))==1
                     data_1(ii,jj,[1 round(size(data_1,3)/2) size(data_1,3)])=[0 0.5 1];
                    end
                end
            end
            data_1=interp1(t_1,permute(data_1,[3 1 2]),t_interp_1,'spline');
        end
    end
    dataF=permute(data_1,[2,3,1]);
    t=t_interp_1;
    clear data_1 t_interp_1
    if median(diff(t))>1e-3
        t=t*1e-3; % make it seconds instead of ms
    end
    % if abs(median(diff(lat)))>1e-3
    %     lat=lat*1e-3; % make it mm instead of meter
    % 
    % end
    % if median(diff(axial))>1e-3
    %     axial=axial*1e-3; % make it mm instead of meter
    % end
    dx=factor.dx;
    % interpolate data
    [xx,yy]=meshgrid(lat,axial);
    disp(axial(1));
    disp(axial(end));
    if isrow(axial)
        axial = axial';
    end
    axialI=axial(1):dx:(axial(end)); % to be safe side

    if median(diff(lat))>0
        latI=lat(1):dx:(lat(end)+.0);
    else
        latI=lat(1):-dx:(lat(end)+.0);
    end
    [xI,yI]=meshgrid(latI,axialI);
    dataInterp=zeros([size(xI),size(dataF,3)]);
    for ii=1:size(dataF,3)
        dataInterp(:,:,ii)=interp2(xx,yy,dataF(:,:,ii),xI,yI,'linear');
    end
    lat=latI'; clear latI
    axial=axialI';clear axialI
    data=dataInterp; clear dataInterp