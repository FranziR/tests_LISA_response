clear; clc; clf
%% Read in FastGB data
[phaseShift, xi, n,dApsi, parameters, tDS]=read_data();
T=max(tDS)+(tDS(2)-tDS(1));
freq0=floor(parameters(:,6).*T)./T;
%% Generate analytic model data
[pr_,ps_,n_,L]=lisa_geometry_modified(tDS);
[YY_, dApsi_, phaseShift_, xi_, ep_, ec_]=lisa_gen_events_modified(tDS,pr_,ps_,n_,L,parameters,freq0);
%% Compare
M=numel(parameters)/8;
freqs=(-numel(tDS)/2:numel(tDS)/2-1)./T;
psFFT=fftshift(fft(phaseShift,[],1),1);
psFFT_=fftshift(fft(phaseShift_,[],1),1);

for jj=1:2
    for ii=1:3
        clf
        err=norm(dApsi_(:,ii,jj)-dApsi(:,ii,jj),'fro')/norm(dApsi_(:,ii,jj),'fro');

        subplot(3,1,1)
        plot(tDS, real(dApsi(:,ii,jj))); hold on
        plot(tDS, real(dApsi_(:,ii,jj)),'--');
        plot(tDS, real(dApsi(:,ii,jj)-dApsi_(:,ii,jj)),'g:');
        title(num2str(err))

        subplot(3,1,2)
        plot(tDS, imag(dApsi(:,ii,jj))); hold on
        plot(tDS, imag(dApsi_(:,ii,jj)),'--');
        plot(tDS, imag(dApsi(:,ii,jj)-dApsi_(:,ii,jj)),'g:');

        subplot(3,1,3)
        plot(tDS, abs(dApsi(:,ii,jj))); hold on
        plot(tDS, abs(dApsi_(:,ii,jj)),'--');
        plot(tDS, abs(dApsi(:,ii,jj)-dApsi_(:,ii,jj)),'g:');

        drawnow
    end
end

for jj=1:3
    nn_=n_(:,:,jj);
    nn=n(:,:,jj);
    for ii=1:3
        clf
        err=norm(nn_(:,ii)-nn(:,ii),'fro')/norm(nn_(:,ii),'fro');

        subplot(3,1,1)
        plot(tDS, real(nn(:,ii))); hold on
        plot(tDS, real(nn_(:,ii)),'--');
        plot(tDS, real(nn(:,ii)-nn_(:,ii)),'g:');
        title(num2str(err))

        subplot(3,1,2)
        plot(tDS, imag(nn(:,ii))); hold on
        plot(tDS, imag(nn_(:,ii)),'--');
        plot(tDS, imag(nn(:,ii)-nn_(:,ii)),'g:');

        subplot(3,1,3)
        plot(tDS, abs(nn(:,ii))); hold on
        plot(tDS, abs(nn_(:,ii)),'--');
        plot(tDS, abs(nn(:,ii)-nn_(:,ii)),'g:');

        drawnow
    end
    drawnow
end

for ii=1:3
    clf
    err=norm(xi_(:,ii)-xi(:,ii),'fro')/norm(xi_(:,ii),'fro');

    subplot(3,1,1)
    plot(tDS, real(xi(:,ii))); hold on
    plot(tDS, real(xi_(:,ii)),'--');
    plot(tDS, real(xi(:,ii)-xi_(:,ii)),'g:');
    title(num2str(err))

    subplot(3,1,2)
    plot(tDS, imag(xi(:,ii))); hold on
    plot(tDS, imag(xi_(:,ii)),'--');
    plot(tDS, imag(xi(:,ii)-xi_(:,ii)),'g:');

    subplot(3,1,3)
    plot(tDS, abs(xi(:,ii))); hold on
    plot(tDS, abs(xi_(:,ii)),'--');
    plot(tDS, abs(xi(:,ii)-xi_(:,ii)),'g:');

    drawnow
end

for jj=1:5
    for ii=1:3
        clf
        err=norm(phaseShift(:,ii,jj)-phaseShift_(:,ii,jj),'fro')/norm(phaseShift(:,ii,jj),'fro');

        subplot(3,2,1)
        plot(tDS, real(phaseShift(:,ii,jj))); hold on
        plot(tDS, real(phaseShift_(:,ii,jj)),'--');
        plot(tDS, real(phaseShift(:,ii,jj)-phaseShift_(:,ii,jj)),'g:');
        title(num2str(err))

        subplot(3,2,3)
        plot(tDS, imag(phaseShift(:,ii,jj))); hold on
        plot(tDS, imag(phaseShift_(:,ii,jj)),'--');
        plot(tDS, imag(phaseShift(:,ii,jj)-phaseShift_(:,ii,jj)),'g:');

        subplot(3,2,5)
        plot(tDS, abs(phaseShift(:,ii,jj))); hold on
        plot(tDS, abs(phaseShift_(:,ii,jj)),'--');
        plot(tDS, abs(phaseShift(:,ii,jj)-phaseShift_(:,ii,jj)),'g:');

        subplot(3,2,2)
        plot(freqs, real(psFFT(:,ii,jj))); hold on
        plot(freqs, real(psFFT_(:,ii,jj)),'--');
        plot(freqs, real(psFFT(:,ii,jj)-psFFT_(:,ii,jj)),'g:');
        xlim([-1e-6 1e-6])
        title(num2str(err))

        subplot(3,2,4)
        plot(freqs, imag(psFFT(:,ii,jj))); hold on
        plot(freqs, imag(psFFT_(:,ii,jj)),'--');
        plot(freqs, imag(psFFT(:,ii,jj)-psFFT_(:,ii,jj)),'g:');
        xlim([-1e-6 1e-6])

        subplot(3,2,6)
        plot(freqs, abs(psFFT(:,ii,jj))); hold on
        plot(freqs, abs(psFFT_(:,ii,jj)),'--');
        plot(freqs, abs(psFFT(:,ii,jj)-psFFT_(:,ii,jj)),'g:');
        xlim([-1e-6 1e-6])

        drawnow
    end
end
%% functions
function [phaseShift, xi, n, dApsi, parameters, t]=read_data()
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

    phaseShift=zeros(numel(t),3,M);
    
    fileID = fopen(fullfile(file_path,'phaseShift12.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    phaseShift(:,1,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);
    
    fileID = fopen(fullfile(file_path,'phaseShift23.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    phaseShift(:,2,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);
    
    fileID = fopen(fullfile(file_path,'phaseShift31.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    phaseShift(:,3,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);

    xi=zeros(numel(t),3,M);   
    fileID = fopen(fullfile(file_path,'xi1.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    xi(:,1,:)=tmp;
    fclose(fileID);
    
    fileID = fopen(fullfile(file_path,'xi2.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    xi(:,2,:)=tmp;
    fclose(fileID);
    
    fileID = fopen(fullfile(file_path,'xi3.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    xi(:,3,:)=tmp;
    fclose(fileID);

    n=zeros(numel(t),3,3);

    fileID = fopen(fullfile(file_path,'n12.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    n(:,1,1)=tmp(1:3:end-2,1);
    n(:,2,1)=tmp(2:3:end-1,1);
    n(:,3,1)=tmp(3:3:end,1);
    fclose(fileID);

    fileID = fopen(fullfile(file_path,'n23.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    n(:,1,2)=tmp(1:3:end-2,1);
    n(:,2,2)=tmp(2:3:end-1,1);
    n(:,3,2)=tmp(3:3:end,1);
    fclose(fileID);

    fileID = fopen(fullfile(file_path,'n31.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    n(:,1,3)=tmp(1:3:end-2,1);
    n(:,2,3)=tmp(2:3:end-1,1);
    n(:,3,3)=tmp(3:3:end,1);
    fclose(fileID);

    dApsi=zeros(numel(t),3,M);
    
    fileID = fopen(fullfile(file_path,'dApsi12.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    dApsi(:,1,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);
    
    fileID = fopen(fullfile(file_path,'dApsi23.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    dApsi(:,2,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);
    
    fileID = fopen(fullfile(file_path,'dApsi31.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    dApsi(:,3,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);

end