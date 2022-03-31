clear; clc; clf
%% Read in FastGB data
[YYFFT, parameters, tDS]=read_data();
T=max(tDS)+(tDS(2)-tDS(1));
freq0=floor(parameters(:,6).*T)./T;
M=numel(parameters)/8;
%% Generate analytic model data
[pr,ps,n,L]=lisa_geometry_modified(tDS);
[~, YYFFT_,E,k,h,fctr]=lisa_gen_events_modified(tDS,pr,ps,n,L,parameters,freq0,T);
%% Generate TLA
ZZFFT_=generate_TLA(freq0,pr,n,L,h,E,k,fctr);
%% Visual comparison of first two signals
freqs=(-numel(tDS)/2:numel(tDS)/2-1)./T;
for jj=1:2
    for ii=1:3
        figure(1)
        clf
        % first column: error reference model (adapted to LDC) - FastGB
        % implementation of LDC
        errAna_FastGB=norm(YYFFT(:,ii)-YYFFT_(:,ii),'fro')/norm(YYFFT(:,ii));
        subplot(2,2,1)
        plot(freqs, real(YYFFT_(:,ii))); hold on
        plot(freqs, real(YYFFT(:,ii)),'--');
        plot(freqs, real(YYFFT_(:,ii)-YYFFT(:,ii)),'g:');
        
        ylabel(['$\Re({X}$(',num2str(ii),'))'], 'Interpreter','latex','FontSize',16)
        xlabel('$\omega$ (Hz)', 'Interpreter','latex','FontSize',16)
        ax=gca; ax.FontSize = 16; ax.TickLabelInterpreter='latex';
        title(['$\mathrm{Reference}_{\mathrm{LDC}}$ vs $\mathrm{FastGB}_{\mathrm{LDC}}$: ', num2str(errAna_FastGB)],'FontSize',16, 'Interpreter','latex')
        grid on
        xlim([-1e-6 1e-6])

        subplot(2,2,3)
        plot(freqs, imag(YYFFT_(:,ii))); hold on
        plot(freqs, imag(YYFFT(:,ii)),'--');
        plot(freqs, imag(YYFFT_(:,ii)-YYFFT(:,ii)),'g:');
        ylabel(['$\Im({X}$(',num2str(ii),'))'], 'Interpreter','latex','FontSize',16)
        xlabel('$\omega$ (Hz)', 'Interpreter','latex','FontSize',16)
        grid on
        xlim([-1e-6 1e-6])
        
        % second column: error reference model (adapted to LDC) - TLA
        % (adapted to LDC)
        errAna_TLA=norm(YYFFT_(:,ii)-ZZFFT_(:,ii),'fro')/norm(YYFFT_(:,ii));
        subplot(2,2,2)
        plot(freqs, real(YYFFT_(:,ii))); hold on
        plot(freqs, real(ZZFFT_(:,ii)),'--');
        plot(freqs, real(YYFFT_(:,ii)-ZZFFT_(:,ii)),'g:');

        ylabel(['$\Re({X}$(',num2str(ii),'))'], 'Interpreter','latex','FontSize',16)
        xlabel('$\omega$ (Hz)', 'Interpreter','latex','FontSize',16)
        ax=gca; ax.FontSize = 16; ax.TickLabelInterpreter='latex';
        title(['$\mathrm{Reference}_{\mathrm{LDC}}$ vs $\mathrm{TLA}_{\mathrm{LDC}}$: ', num2str(errAna_TLA)],'FontSize',16, 'Interpreter','latex')
        grid on
        
        xlim([-1e-6 1e-6])

        subplot(2,2,4)
        plot(freqs, imag(YYFFT_(:,ii))); hold on
        plot(freqs, imag(ZZFFT_(:,ii)),'--');
        plot(freqs, imag(YYFFT_(:,ii)-ZZFFT_(:,ii)),'g:');
        ylabel(['$\Im({X}$(',num2str(ii),'))'], 'Interpreter','latex','FontSize',16)
        xlabel('$\omega$ (Hz)', 'Interpreter','latex','FontSize',16)
        ax=gca; ax.FontSize = 16; ax.TickLabelInterpreter='latex';
        grid on

        xlim([-1e-6 1e-6])
        
        drawnow
    end
end
%% Overview
errorRef_FGB=zeros(M,1);
errorRef_TLA=zeros(M,1);
for mm=1:M
    errorRef_FGB(mm,1)=norm(YYFFT_(:,:,mm)-YYFFT(:,:,mm),'fro')/norm(YYFFT_(:,:,mm),'fro');
    errorRef_TLA(mm,1)=norm(YYFFT_(:,:,mm)-ZZFFT_(:,:,mm),'fro')/norm(YYFFT_(:,:,mm),'fro');
end
figure(2)
clf
semilogy((1:M), errorRef_FGB,'b*','DisplayName','$\mathrm{Reference}_{\mathrm{LDC}}$ vs $\mathrm{FastGB}_{\mathrm{LDC}}$')
hold on
semilogy((1:M), errorRef_TLA,'g*','DisplayName','$\mathrm{Reference}_{\mathrm{LDC}}$ vs $\mathrm{TLA}_{\mathrm{LDC}}$')
ylabel('relative error in $\mathrm{L}_2$ norm', 'Interpreter','latex','FontSize',16)
ax=gca; ax.FontSize = 16; ax.TickLabelInterpreter='latex';
lg=legend(); lg.Interpreter='latex'; lg.Position=[0.7 0.15 0.05 0.05]; lg.NumColumns=2;
title(['relative error in all ', num2str(M),' signals'],'FontSize',16, 'Interpreter','latex')
ylim([1e-10 1e-3])
grid on
%% functions
function ZZFFT_=generate_TLA(freq0, pr, n, L, h, E, k, fctr)
    M=numel(freq0);
    N=size(pr,1);
    ZZ_=zeros(N,3,M);
    for m=1:M
        A=get_A_LDC(freq0(m),pr,n,L);
        tmp=squeeze(sum(sum( h(:,m).*A.*E(:,m).'.*reshape([1;k(:,m)],[1 1 4]),2),3));
        ZZ_(:,:,m)=ifft(ifftshift(fftshift(fft(tmp,[],1),1).*fctr(:,:,m),1),[],1)./(2*N);
        ZZ_(:,:,m)=ZZ_(:,:,m).*exp(pi/2*1i-2*1i*pi*freq0(m).*L(:,[1 2 3])./299792458);
    end
    ZZFFT_=fftshift(fft(ZZ_,[],1),1);
end


function [YY, parameters, t]=read_data()
    Tobs=33554432; N=1024; dt=Tobs/N; t=(0:dt:Tobs-dt)';

    file_path = fileparts(mfilename('fullpath'));

    fileID = fopen(fullfile(file_path,'parameters.bin')); 
    tmp=fread(fileID,'double');
    parameters=permute(reshape(tmp, 8, []),[2,1]);
    fclose(fileID);
    M=numel(parameters)/8;

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

end