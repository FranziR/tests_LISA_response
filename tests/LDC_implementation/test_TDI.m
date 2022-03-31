clear; clc; clf
%% Read in FastGB data
[YY, YYLDC]=read_data();
M=25;
%% Compare
for jj=1:2
    for ii=1:3
        clf
        figure(10)
        err=norm(YYLDC(:,ii,jj)-YY(:,ii,jj),'fro')/norm(YYLDC(:,ii,jj),'fro');

        subplot(3,1,1)
        plot(real(YY(:,ii,jj))); hold on
        plot(real(YYLDC(:,ii,jj)),'--');
        plot(real(YY(:,ii,jj)-YYLDC(:,ii,jj)),'g:');
        title(num2str(err))

        subplot(3,1,2)
        plot(imag(YY(:,ii,jj))); hold on
        plot(imag(YYLDC(:,ii,jj)),'--');
        plot(imag(YY(:,ii,jj)-YYLDC(:,ii,jj)),'g:');

        subplot(3,1,3)
        plot(abs(YY(:,ii,jj))); hold on
        plot(abs(YYLDC(:,ii,jj)),'--');
        plot(abs(YY(:,ii,jj)-YYLDC(:,ii,jj)),'g:');

        drawnow
    end
end

error=zeros(M,1);
for mm=1:M
    error(mm,1)=norm(YYLDC(:,:,mm)-YY(:,:,mm),'fro')/norm(YYLDC(:,:,mm),'fro');
end
figure(2)
clf
semilogy((1:M), error,'b*')
drawnow
%% functions
function [YY, YYLDC]=read_data()
    M=25; N=32768;
    file_path = fileparts(mfilename('fullpath')); 
    
    YY=zeros(N,3,M);
    
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

    YYLDC=zeros(N,3,M);
    
    fileID = fopen(fullfile(file_path,'X_LDC.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    YYLDC(:,1,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);   
    
    fileID = fopen(fullfile(file_path,'Y_LDC.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    YYLDC(:,2,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);
    
    fileID = fopen(fullfile(file_path,'Z_LDC.bin')); 
    tmp=fread(fileID,'double');
    tmp=reshape(tmp, [], M);
    YYLDC(:,3,:)=tmp(1:2:end-1,:)+1i*tmp(2:2:end,:);
    fclose(fileID);

end