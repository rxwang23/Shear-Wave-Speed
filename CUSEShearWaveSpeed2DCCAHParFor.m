function [cs,cc,csx,csy] = CUSEShearWaveSpeed2DCCAHParFor(data,prf,dx,dy,direc,stepsize,winsize,IntFactor)
    hsize = round(winsize/2);
    nump = winsize - stepsize + 1; % number of points to average by AH
    offset = (nump-1)/2; % offset to get the correct index;
    
    %% Turn on parallel processing
    %% Transpose the data if the waves are ud or du
    if strcmp(direc,'ud') || strcmp(direc,'du')
        temp = permute(data,[2 1 3]);
        clear data;
        data = temp;clear temp;
    end
    %% Initiate parameters
    [Nax,Nlat,Nt] = size(data);
    %% interp data
    t_interp=1:(1/IntFactor):Nt;
    data=permute(interp1(1:Nt,permute(data,[3 1 2]),t_interp,'pchip'),[2 3 1]);
    % tukey window
    tWin=tukeywin(size(data,3),.25);
    tWinMat=permute(repmat(tWin,[1 size(data(:,:,1))]),[2 3 1]);
    data=data.*tWinMat;
    % x-direction speed calculation (parallel to wave propagation direction)
    csxtemp = zeros(Nax,Nlat-stepsize);
    ccxtemp = zeros(Nax,Nlat-stepsize);
    csxidx = zeros(Nax,Nlat-stepsize);
    % z-direction speed calculation
    csytemp = zeros(Nax-stepsize,Nlat);
    ccytemp = zeros(Nax-stepsize,Nlat);
    csyidx = zeros(Nax-stepsize,Nlat);
    for i = 1:Nax
        mv = 0;
        idx = 0;
        for j = 1:Nlat
            dataX=squeeze(data(i,j,:));
            % calculate SWV in x direction
            if j<(Nlat-stepsize+1)
                dataX1=squeeze(data(i,j+stepsize,:));
                % Cross-correlation
                [r,lags] = xcorr(dataX,dataX1,'coeff');
                % Find the correct correlation coefficient peak
                if strcmp(direc,'lr') || strcmp(direc,'ud')
                    [mv,idx] = max(r(lags<0));
                end
                if strcmp(direc,'rl') || strcmp(direc,'du')
                    [mv,idx] = max(r(lags>0));
                    idx = idx + (length(lags) - 1)/2 + 1;
                end
                % Calculate time of flight
                del = abs(lags(idx))*(1/(prf*IntFactor));
                % Calculate shear wave speed
                csxtemp(i,j) = stepsize*dx/del;
                ccxtemp(i,j) = mv;
                csxidx(i,j) = (j + j + stepsize)/2;
            end
             if i<(Nax-stepsize+1)
                dataZ1=squeeze(data(i+stepsize,j,:));
                % Cross-correlation
                [r2,lags2] = xcorr(dataX,dataZ1,'coeff');
                
                % Find the correct correlation coefficient peak
                [mv6,idx2] = max(r2);
        %         if (lags2(idx2)== 0)
                if abs(lags2(idx2))<2*IntFactor
    %                 del2 = eps; % if there is no Vy component, make delay an infinitesimal number
    %                 mv6 = 0;
                    csytemp(i,j) = 35;
                    ccytemp(i,j) = mv6;
                else
                    [mv6,idx6] = max(r2);
                    % Calculate time of flight
                    del2 = abs(lags2(idx6))*(1/(prf*IntFactor));
                    csytemp(i,j) = stepsize*dy/del2;
                    ccytemp(i,j) = mv6;
                end
                % Calculate shear wave speed
               
                csyidx(i,j) = (i + i + stepsize)/2;
             end
        end
    end
    
    % csytemp(ccytemp<0.5)=0;
    % csxtemp(ccxtemp<0.5)=0;
    % 
    % disp('Number of columns in Y:');
    % disp(size(Y, 2));
    
    %% Combine the Csx and Csy calculations using Anderssen Hegland's method
    % Form a distance-weighting matrix
    [X,Y] = meshgrid(-offset:offset,-hsize:hsize);
    vxweight = (1./sqrt(X.^2 + Y.^2));
    vxweight(isinf(vxweight)) = max(vxweight(~isinf(vxweight)));
    [X,Y] = meshgrid(-hsize:hsize,-offset:offset);
    vyweight = (1./sqrt(X.^2 + Y.^2));
    vyweight(isinf(vyweight)) = max(vyweight(~isinf(vyweight)));
    csx = zeros(Nax,Nlat);
    csy = zeros(Nax,Nlat);
    cs = zeros(Nax,Nlat);
    cc = zeros(Nax,Nlat);
    % Start the calculation
    for i = 2:Nax-1
        hsize1 = 0;
        mv3 = 0;
        idx3 = 0;
        mv4 = 0;
        idx4 = 0;
        for j = 2:Nlat-1
            if (j>hsize && j<=Nlat-hsize && i>hsize && i<=Nax-hsize)
                xidx = j-offset:j+offset;
                yidx = i-offset:i+offset;
                indx = find(csxidx(i,:) >= xidx(1) & csxidx(i,:) <= xidx(end));
                indy = find(csyidx(:,j) >= yidx(1) & csyidx(:,j) <= yidx(end));
                csxtemp1 = csxtemp(i-hsize:i+hsize,indx);
                ccx = ccxtemp(i-hsize:i+hsize,indx);
                csytemp1 = csytemp(indy,j-hsize:j+hsize);
                ccy = ccytemp(indy,j-hsize:j+hsize);
                csx(i,j) = sum(sum(csxtemp1.*((ccx.^2.*vxweight)./sum(sum(ccx.^2.*vxweight)))));
                csy(i,j) = sum(sum(csytemp1.*((ccy.^2.*vyweight)./sum(sum(ccy.^2.*vyweight)))));
                if isnan(csy(i,j))
                    cs(i,j) = csx(i,j);
                else
                    if nanmedian(ccx(:))>0.5 && nanmedian(ccy(:))>0.5 
                        cs(i,j) = csx(i,j)*csy(i,j)./(sqrt(csx(i,j)^2+csy(i,j)^2));
                        cc(i,j) = min(mean2(ccx),mean2(ccy));
                    elseif nanmedian(ccx(:))>0.5 && nanmedian(ccy(:))<0.5 
                        cs(i,j) = csx(i,j);
                        cc(i,j) =mean2(ccx);
                    elseif nanmedian(ccx(:))<0.5 && nanmedian(ccy(:))>0.5 
                        cs(i,j) = csy(i,j);
                        cc(i,j) =mean2(ccy);
                    else
    %             if nanmedian(ccTempx(:))< nanmedian(ccTempy(:))
    %                 csV(ii,jj)=csy(ii,jj);
    %                 r2Avg(ii,jj)=nanmedian(ccTempy(:));
    %             else
    %                 csV(ii,jj)=csx(ii,jj);
    %                 r2Avg(ii,jj)=nanmedian(ccTempx(:));
    %             end
                    
                    end
    %                 cs(i,j) = csx(i,j)*csy(i,j)./(sqrt(csx(i,j)^2+csy(i,j)^2));
                end
    %             cc(i,j) = min(nanmedian(ccx(:)),nanmedian(ccy(:)));
    %             cc(i,j) = min(mean2(ccx),mean2(ccy));
            end
            if (i<=hsize || i>Nax-hsize)
                if (j<=hsize)
                    hsize1 = j-1;
                else if (j>Nlat-hsize)
                        hsize1 = Nlat-j;
                else if (j>hsize && j<=Nlat-hsize)
                        hsize1 = hsize;
                    end
                    end
                end
                temp5 = squeeze(data(i,j-hsize1,:));
                temp6 = squeeze(data(i,j+hsize1,:));
                temp5 = interp(temp5,IntFactor);
                temp6 = interp(temp6,IntFactor);
                win3 = tukeywin(length(temp5),.25);
                temp5 = temp5.*win3;
                temp6 = temp6.*win3;
                % Cross-correlation
                [r3,lags3] = xcorr(temp5,temp6,'coeff');
                % Find the correct correlation coefficient peak
                if strcmp(direc,'lr') || strcmp(direc,'ud') 
                    [mv3,idx3] = max(r3(lags3<0));
                end
                if strcmp(direc,'rl') || strcmp(direc,'du') 
                    [mv3,idx3] = max(r3(lags3>0));
                    idx3 = idx3 + (length(lags3) - 1)/2 + 1;
                end
                % Calculate time of flight
                del3 = abs(lags3(idx3))*(1/(prf*IntFactor));
                cs(i,j) = (hsize1*2)*dx/del3;
                cc(i,j) = mv3;
                csx(i,j) = cs(i,j);
                csy(i,j) = 0;
            end
            if (i>hsize && i<=Nax-hsize && (j<=hsize || j>Nlat-hsize))
                if (j<=hsize)
                    hsize1 = j-1;
                else if (j>Nlat-hsize)
                        hsize1 = Nlat-j;
                    end
                end
                temp7 = squeeze(data(i,j-hsize1,:));
                temp8 = squeeze(data(i,j+hsize1,:));
                temp7 = interp(temp7,IntFactor);
                temp8 = interp(temp8,IntFactor);
                win4 = tukeywin(length(temp7),.25);
                temp7 = temp7.*win4;
                temp8 = temp8.*win4;
                % Cross-correlation
                [r4,lags4] = xcorr(temp7,temp8,'coeff');
                % Find the correct correlation coefficient peak
                if strcmp(direc,'lr') || strcmp(direc,'ud')
                    [mv4,idx4] = max(r4(lags4<0));
                end
                if strcmp(direc,'rl') || strcmp(direc,'du')
                    [mv4,idx4] = max(r4(lags4>0));
                    idx4 = idx4 + (length(lags4) - 1)/2 + 1;
                end
                % Calculate time of flight
                del4 = abs(lags4(idx4))*(1/(prf*IntFactor));
                cs(i,j) = (hsize1*2)*dx/del4;
                cc(i,j) = mv4;
                csx(i,j) = cs(i,j);
                csy(i,j) = 0;
            end
        end
    end
    %% Pack the edge data to original dimension
    csx(:,1) = csx(:,2);
    csx(:,Nlat) = csx(:,Nlat-1);
    csx(1,:) = csx(2,:);
    csx(Nax,:) = csx(Nax-1,:);
    
    csy(:,1) = csy(:,2);
    csy(:,Nlat) = csy(:,Nlat-1);
    csy(1,:) = csy(2,:);
    csy(Nax,:) = csy(Nax-1,:);
    
    cs(:,1) = cs(:,2);
    cs(:,Nlat) = cs(:,Nlat-1);
    cs(1,:) = cs(2,:);
    cs(Nax,:) = cs(Nax-1,:);
    
    cc(:,1) = cc(:,2);
    cc(:,Nlat) = cc(:,Nlat-1);
    cc(1,:) = cc(2,:);
    cc(Nax,:) = cc(Nax-1,:);
    %%
    if strcmp(direc,'ud') || strcmp(direc,'du')
        cs = cs';
        cc = cc';
    end
