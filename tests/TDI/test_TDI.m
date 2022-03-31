clear; clc; clf
%% Read in FastGB data
[YY, parameters, tDS]=read_data();
T=max(tDS)+(tDS(2)-tDS(1));
freq0=floor(parameters(:,6).*T)./T;
%% Generate analytic model data
[pr_,ps_,n_,L]=lisa_geometry_modified(tDS);
[~, YY_,E,k,h,fctr]=lisa_gen_events_modified(tDS,pr_,ps_,n_,L,parameters,freq0,T);
%% Compare
M=numel(parameters)/8;
freqs=(-numel(tDS)/2:numel(tDS)/2-1)./T;

for jj=1:M
    for ii=1:3
        clf
        err=norm(YY_(:,ii,jj)-YY(:,ii,jj),'fro')/norm(YY_(:,ii,jj),'fro');

        subplot(3,1,1)
        plot(freqs, real(YY(:,ii,jj))); hold on
        plot(freqs, real(YY_(:,ii,jj)),'--');
        plot(freqs, real(YY(:,ii,jj)-YY_(:,ii,jj)),'g:');
        title(num2str(err))

        subplot(3,1,2)
        plot(freqs, imag(YY(:,ii,jj))); hold on
        plot(freqs, imag(YY_(:,ii,jj)),'--');
        plot(freqs, imag(YY(:,ii,jj)-YY_(:,ii,jj)),'g:');

        subplot(3,1,3)
        plot(freqs, abs(YY(:,ii,jj))); hold on
        plot(freqs, abs(YY_(:,ii,jj)),'--');
        plot(freqs, abs(YY(:,ii,jj)-YY_(:,ii,jj)),'g:');

        drawnow
    end
end
%% functions
function [YY, parameters, t]=read_data()
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