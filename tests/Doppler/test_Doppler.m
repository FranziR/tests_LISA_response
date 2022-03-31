clear; clc; clf
%% Read in FastGB data
[YY, YYFFT, parameters, tDS]=read_data();
T=max(tDS)+(tDS(2)-tDS(1));
freq0=floor(parameters(:,6).*T)./T;
%% Generate analytic model data
[pr_,ps_,n_,L]=lisa_geometry_modified(tDS);
[YY_]=lisa_gen_events_modified(tDS,ps_,n_,L,parameters,freq0);
%% Compare
M=numel(parameters)/8;
freqs=(-numel(tDS)/2:numel(tDS)/2-1)./T;
YFFT_=fftshift(fft(YY_,[],1),1);

for jj=1:M
    for ii=1:6
        clf
        err=norm(YY_(:,ii,jj)-YY(:,ii,jj),'fro')/norm(YY_(:,ii,jj),'fro');

        subplot(3,1,1)
        plot(tDS, real(YY(:,ii,jj))); hold on
        plot(tDS, real(YY_(:,ii,jj)),'--');
        plot(tDS, real(YY(:,ii,jj)-YY_(:,ii,jj)),'g:');
        title(num2str(err))

        subplot(3,1,2)
        plot(tDS, imag(YY(:,ii,jj))); hold on
        plot(tDS, imag(YY_(:,ii,jj)),'--');
        plot(tDS, imag(YY(:,ii,jj)-YY_(:,ii,jj)),'g:');

        subplot(3,1,3)
        plot(tDS, abs(YY(:,ii,jj))); hold on
        plot(tDS, abs(YY_(:,ii,jj)),'--');
        plot(tDS, abs(YY(:,ii,jj)-YY_(:,ii,jj)),'g:');

        drawnow
    end
end

for jj=1:M
    for ii=1:6
        clf
        err=norm(YFFT_(:,ii,jj)-YYFFT(:,ii,jj),'fro')/norm(YFFT_(:,ii,jj),'fro');

        subplot(3,1,1)
        plot(freqs, real(YYFFT(:,ii,jj))); hold on
        plot(freqs, real(YFFT_(:,ii,jj)),'--');
        plot(freqs, real(YYFFT(:,ii,jj)-YFFT_(:,ii,jj)),'g:');
        title(num2str(err))

        subplot(3,1,2)
        plot(freqs, imag(YYFFT(:,ii,jj))); hold on
        plot(freqs, imag(YFFT_(:,ii,jj)),'--');
        plot(freqs, imag(YYFFT(:,ii,jj)-YFFT_(:,ii,jj)),'g:');

        subplot(3,1,3)
        plot(freqs, abs(YYFFT(:,ii,jj))); hold on
        plot(freqs, abs(YFFT_(:,ii,jj)),'--');
        plot(freqs, abs(YYFFT(:,ii,jj)-YFFT_(:,ii,jj)),'g:');

        drawnow
    end
end
%% functions
function [y, yFFT, parameters, t]=read_data()
    Tobs=33554432;
    N=1024; 
    dt=Tobs/N;
    t=(0:dt:Tobs-dt)';

    file_path = fileparts(mfilename('fullpath')); 

    fileID = fopen(fullfile(file_path,'parameters.bin')); 
    tmp=fread(fileID,'double');
    parameters=permute(reshape(tmp, 8, []),[2,1]);
    fclose(fileID);
    M=numel(parameters)/8;

    y=zeros(numel(t),6,M);
    
    fileID = fopen(fullfile(file_path,'y12.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    y(:,1,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);
    
    fileID = fopen(fullfile(file_path,'y23.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    y(:,2,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);
    
    fileID = fopen(fullfile(file_path,'y31.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    y(:,3,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);

    fileID = fopen(fullfile(file_path,'y21.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    y(:,4,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);

    fileID = fopen(fullfile(file_path,'y32.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    y(:,5,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);

    fileID = fopen(fullfile(file_path,'y13.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    y(:,6,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);

    yFFT=zeros(numel(t),6,M);
    
    fileID = fopen(fullfile(file_path,'yFFT12.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    yFFT(:,1,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);
    
    fileID = fopen(fullfile(file_path,'yFFT23.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    yFFT(:,2,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);
    
    fileID = fopen(fullfile(file_path,'yFFT31.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    yFFT(:,3,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);

    fileID = fopen(fullfile(file_path,'yFFT21.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    yFFT(:,4,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);

    fileID = fopen(fullfile(file_path,'yFFT32.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    yFFT(:,5,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);

    fileID = fopen(fullfile(file_path,'yFFT13.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    yFFT(:,6,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);

    yFFT=fftshift(yFFT,1);

end