clear; clc; clf
%% Read in FastGB data
[~, YYFFT, parameters, tDS]=read_data();
T=max(tDS)+(tDS(2)-tDS(1));
freq0=floor(parameters(6).*T)./T;
M=1;
%% Generate analytic model data
[pr,ps,n,L]=lisa_geometry_modified(tDS);
[~, YYFFT_,E,k,h,fctr]=lisa_gen_events_modified(tDS,pr,ps,n,L,parameters,freq0,T);
YY_=ifft(ifftshift(YYFFT_,1),[],1);
%% Generate TLA
ZZ_=zeros(length(tDS),3,M);
for m=1:M
    A=get_A_LDC(freq0(m),pr,n,L);
    tmp=squeeze(sum(sum( h(:,m).*A.*E(:,m).'.*reshape([1;k(:,m)],[1 1 4]),2),3));
    ZZ_(:,:,m)=ifft(ifftshift(fftshift(fft(tmp,[],1),1).*fctr(:,:,m),1),[],1)./(2*numel(tDS));
    ZZ_(:,:,m)=ZZ_(:,:,m).*exp(pi/2*1i-2*1i*pi*parameters(m,6).*L(:,[1 2 3])./299792458);
end
ZZFFT_=fftshift(fft(ZZ_,[],1),1);
%% Compare
freqs=(-numel(tDS)/2:numel(tDS)/2-1)./T;
for jj=1:1
    for ii=1:3
        clf
        % first column: error analytic - TLA
%         errAna_TLA=norm(YYFFT_(:,ii,jj)-ZZFFT_(:,ii,jj),'fro')/norm(YYFFT_(:,ii,jj));
        errAna_TLA=norm(YYFFT_(:,ii)-ZZFFT_(:,ii),'fro')/norm(YYFFT_(:,ii));
        subplot(3,2,1)
        plot(freqs, real(YYFFT_(:,ii))); hold on
        plot(freqs, real(ZZFFT_(:,ii)),'--');
        plot(freqs, real(YYFFT_(:,ii)-ZZFFT_(:,ii)),'g:');
        xlim([-1e-6 1e-6])
        title(['Analytic vs TLA: ', num2str(errAna_TLA)])

        subplot(3,2,3)
        plot(freqs, imag(YYFFT_(:,ii))); hold on
        plot(freqs, imag(ZZFFT_(:,ii)),'--');
        plot(freqs, imag(YYFFT_(:,ii)-ZZFFT_(:,ii)),'g:');
        xlim([-1e-6 1e-6])

        subplot(3,2,5)
        plot(freqs, abs(YYFFT_(:,ii))); hold on
        plot(freqs, abs(ZZFFT_(:,ii)),'--');
        plot(freqs, abs(YYFFT_(:,ii)-ZZFFT_(:,ii)),'g:');
        xlim([-1e-6 1e-6])
        
        % second column: error analytic - FastGB
        errAna_FastGB=norm(YYFFT(:,ii)-YYFFT_(:,ii),'fro')/norm(YYFFT(:,ii));
        subplot(3,2,2)
        plot(freqs, real(YYFFT_(:,ii))); hold on
        plot(freqs, real(YYFFT(:,ii)),'--');
        plot(freqs, real(YYFFT_(:,ii)-YYFFT(:,ii)),'g:');
        xlim([-1e-6 1e-6])
        title(['Analytic vs FastGB: ', num2str(errAna_FastGB)])

        subplot(3,2,4)
        plot(freqs, imag(YYFFT_(:,ii))); hold on
        plot(freqs, imag(YYFFT(:,ii)),'--');
        plot(freqs, imag(YYFFT_(:,ii)-YYFFT(:,ii)),'g:');
        xlim([-1e-6 1e-6])

        subplot(3,2,6)
        plot(freqs, abs(YYFFT_(:,ii))); hold on
        plot(freqs, abs(YYFFT(:,ii)),'--');
        plot(freqs, abs(YYFFT_(:,ii)-YYFFT(:,ii)),'g:');
        xlim([-1e-6 1e-6])

        drawnow
    end
end
%% functions
function [YY, YYSL, parameters, t]=read_data()
    Tobs=33554432;
    N=1024; 
    dt=Tobs/N;
    t=(0:dt:Tobs-dt)';

    file_path = fileparts(mfilename('fullpath'));
    M=1; 

    fileID = fopen('parameters.txt','r');
    tmp=fscanf(fileID,'%f %f %f %f %f %f %f %f');
    fclose(fileID);
    parameters=permute(reshape(tmp, 8, 10),[2, 1]);
    parameters=parameters(1,:);

    YY=zeros(numel(t),3,M);
    
    fileID = fopen(fullfile(file_path,'X.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    YY(:,1,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);   
    
    fileID = fopen(fullfile(file_path,'Y.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    YY(:,2,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);
    
    fileID = fopen(fullfile(file_path,'Z.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    YY(:,3,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);
    
    YYSL=zeros(numel(t),3,M);
    
    fileID = fopen(fullfile(file_path,'XSL.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    YYSL(:,1,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);   
    
    fileID = fopen(fullfile(file_path,'YSL.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    YYSL(:,2,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);
    
    fileID = fopen(fullfile(file_path,'ZSL.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    YYSL(:,3,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);

end