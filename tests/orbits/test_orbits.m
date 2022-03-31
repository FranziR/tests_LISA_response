clear all; clc; clf;
%% 0. read in FastGB data
[Tobs, t, p]=reading_data(); 
%% 1. generate original data
[~,p_,~,~]=lisa_geometry_modified(t);

clf
for jj=1:3
    subplot(1,3,jj)
    plot(t, p_(:,jj,2)); hold on
    plot(t, p(:,jj),'--');
    drawnow
end
error=norm(p_(:,:,2)-p,'fro')/norm(p_(:,:,2),'fro')
%% functions
function [Tobs, t, p]=reading_data()
    Tobs=33554432;
    N=1024; 
    dt=Tobs/N;
    t=(0:dt:Tobs-dt);

    p=zeros(numel(t),3);
    % read time
    filepath = fileparts(mfilename('fullpath'));
    fileID=fopen(fullfile(filepath,'p1_x.bin')); 
    p(:,1)=fread(fileID,'double'); fclose(fileID);

    fileID=fopen(fullfile(filepath,'p1_y.bin')); 
    p(:,2)=fread(fileID,'double'); fclose(fileID);

    fileID=fopen(fullfile(filepath,'p1_z.bin')); 
    p(:,3)=fread(fileID,'double'); fclose(fileID);
    
end

function data0=rifft_and_interpolation(f, T, N, dataFFT)
    ND=numel(dataFFT)/3;
    kmin=round(f*T-ND/2);
    X=zeros(N,3);
    if f*T-round(f*T)>0
        X(kmin+N/2+1:kmin+ND+N/2,:) = dataFFT; % dataset 1,3,9,10
        X(N/2-kmin-ND:N/2-kmin-1,:) = conj(flip(dataFFT,1));
    else
        X(kmin+N/2:kmin-1+ND+N/2,:) = dataFFT; % dataset 2,4,5,6,7,8
        X(N/2-kmin-ND+1:N/2-kmin,:) = conj(flip(dataFFT,1));
    end
    data0=(N/ND)*ifft(ifftshift(X,1),[],1);
end

function data=downsampling(data0, t, R, freq0)
    data2=ifft(fft(data0).*[1;2*ones(size(data0,1)/2-1,1);1;zeros(size(data0,1)/2-1,1)]);
    data3=data2.*exp(-2*pi*1i*t.*reshape(freq0.*ones(3,1),[1,3]));
    data=1/R*data3(1:R:end,:);
end
