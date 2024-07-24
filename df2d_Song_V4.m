function [dfiltered, dfilters] = df2d_Song_V4(data, dx, dy, cutoffs, order, power, angles)
% [dfiltered,dfilters] = df2d_Song_V4(data, dx, dy, cutoffs, order, power, varargin);
% performs 2D directional filtering in desired directions
% This version only works for the data with same axial and lateral resolutions

% data - input data
% dx - lateral resolution (m)
% dy - axial resolution (m)
% cutoffs - cutoffs for butterworth bp filter, normalized to 1/2 of spatial
            % sampling frequency (e.g., [.1 .5] = [.1*1/resol/2 .5*1/resol/2] wavenumber
% order - order of butterworth bp filter (actual order = 2*input order)
% power - power of dot product, must be even, higher powers have sharper directional bw
% angles - input angles of the desired filtering directions, in format of [theta1, theta2, ..., thetaN],theta in radians

%dfiltered - returns filtered data in a 4-D matrix, the 4th dimension corresponds to input angles
%dfilters: returns the 2D plots of the directional filters with angles specified by the input angles
%%
warning('off','MATLAB:log:logOfZero');
%% Load info of desired directions
nofilts = length(angles);

%% windowing the data to remove ringing
nx=size(data,2);ny=size(data,1); nt = size(data,3);

%% Center the data
dcShiftx = 1 - mod(nx,2)/2;
dcShifty = 1 - mod(ny,2)/2;

%% load the high and low cutoff wavenumbers
hc=cutoffs(2)*(1/dx/2);
lc=cutoffs(1)*(1/dx/2);

%% Define a grid to calculate cosine angles
[ii,jj] = ndgrid(single(1:ny),single(1:nx));
rhoVector = cat(3, (1/dy/ny)*(ii-(ny/2+dcShifty)), (1/dx/nx)*(jj-(nx/2+dcShiftx)));
% Get the distance matrix
rho=sqrt(sum(rhoVector.^2,3));
% Define bandpass filter
bpFilter=(1+(((hc*lc)-rho.^2)./((hc-lc)*rho)).^(2*order)).^-.5;

%% Calculate the filters
% For positive temporal frequencies
dfVector = zeros(2,nofilts);
for n = 1:nofilts
    dfVector(:,n) = [-sin(angles(n)),cos(angles(n))];
end

% For negative temporal frequencies
dfVector2 = zeros(2,nofilts);
for n = 1:nofilts
    dfVector2(:,n) = [sin(angles(n)),-cos(angles(n))];
end

normRhoVector = rhoVector./repmat(rho,[1 1 2]);
% set dc to 0
ind= isnan(normRhoVector); normRhoVector(ind)=0;
% Get the filter (df: for positive frequencies; dfc: for negative
% frequencies)
for dirn = 1:nofilts
    temp = dfVector(1,dirn)*normRhoVector(:,:,1) + dfVector(2,dirn)*normRhoVector(:,:,2);
    ind= temp<0; temp(ind)=0;
    df(:,:,dirn)=(temp.^power).*bpFilter;
end

for dirn = 1:nofilts
    temp = dfVector2(1,dirn)*normRhoVector(:,:,1) + dfVector2(2,dirn)*normRhoVector(:,:,2);
    ind=find(temp<0); temp(ind)=0;
    dfc(:,:,dirn)=(temp.^power).*bpFilter;
end
     
dfilters=df;

%% Apply the filter
xf=fft(data,[],3); 
% find the checkpoint where the directional filter direction is rotated by
% 180 degrees
if rem(nt,2)
    checkpoint = (nt+1)/2;
else
    checkpoint = nt/2;
end
% apply the filter
nf = size(xf,3);
kspdf = zeros(ny,nx,nt,nofilts);
for n = 1:nf
    if n <= checkpoint
        ksp = fftshift(fft2(xf(:,:,n)));
        kspdf(:,:,n,:) = reshape(repmat(ksp,[1 1 nofilts]).*dfc,ny,nx,1,nofilts);
    else
        ksp = fftshift(fft2(xf(:,:,n)));
        kspdf(:,:,n,:) = reshape(repmat(ksp,[1 1 nofilts]).*df,ny,nx,1,nofilts);
    end
end

%% Return the filtered the data
dfiltered = real(ifft(ifft(ifft(ifftshift(ifftshift(kspdf,1),2),[],1),[],2),[],3));

warning('on','MATLAB:log:logOfZero');
