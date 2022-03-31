clear; clc; clf
%% Read in FastGB data
[YY, YYSL, parameters, tDS]=read_data();
T=max(tDS)+(tDS(2)-tDS(1));
freq0=floor(parameters(6).*T)./T;
%% Generate analytic model data
[pr_,ps_,n_,L]=lisa_geometry_modified(tDS);
[YY_, YYSL_,E,k,h,fctr]=lisa_gen_events_modified(tDS,pr_,ps_,n_,L,parameters,freq0,T);
% YY_=fftshift(fft(YY_,[],1),1);
% YYSL_=fftshift(fft(YYSL_,[],1),1);
%% Compare
M=numel(parameters)/8;
freqs=(-numel(tDS)/2:numel(tDS)/2-1)./T;

for ii=1:3
    clf
    err=norm(YY_(:,ii)-YY(:,ii),'fro')/norm(YY_(:,ii),'fro');

    subplot(3,1,1)
    plot(freqs, real(YY(:,ii))); hold on
    plot(freqs, real(YY_(:,ii)),'--');
    plot(freqs, real(YY(:,ii)-YY_(:,ii)),'g:');
    title(num2str(err))

    subplot(3,1,2)
    plot(freqs, imag(YY(:,ii))); hold on
    plot(freqs, imag(YY_(:,ii)),'--');
    plot(freqs, imag(YY(:,ii)-YY_(:,ii)),'g:');

    subplot(3,1,3)
    plot(freqs, abs(YY(:,ii))); hold on
    plot(freqs, abs(YY_(:,ii)),'--');
    plot(freqs, abs(YY(:,ii)-YY_(:,ii)),'g:');

    drawnow
end


for ii=1:3
    clf
    err=norm(YYSL_(:,ii)-YYSL(:,ii),'fro')/norm(YYSL_(:,ii),'fro');

    subplot(3,1,1)
    plot(freqs, real(YYSL(:,ii))); hold on
    plot(freqs, real(YYSL_(:,ii)),'--');
    plot(freqs, real(YYSL(:,ii)-YYSL_(:,ii)),'g:');
    title(num2str(err))

    subplot(3,1,2)
    plot(freqs, imag(YYSL(:,ii))); hold on
    plot(freqs, imag(YYSL_(:,ii)),'--');
    plot(freqs, imag(YYSL(:,ii)-YYSL_(:,ii)),'g:');

    subplot(3,1,3)
    plot(freqs, abs(YYSL(:,ii))); hold on
    plot(freqs, abs(YYSL_(:,ii)),'--');
    plot(freqs, abs(YYSL(:,ii)-YYSL_(:,ii)),'g:');

    drawnow
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