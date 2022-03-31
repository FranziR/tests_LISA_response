function A=get_A_LDC(freq,p,n,L)
    p=p./299792458;
    L=L./299792458;
    Nt=size(L,1);
    nE=zeros(Nt,9,3);nE_tot=zeros(Nt,6,6);
    for j=1:3 
        tmp=reshape(permute(n(:,:,j),[2 1 3]),[3 1 Nt]);
        nE(:,:,j)=reshape(tmp.*permute(tmp,[2 1 3]),[9 Nt])';
    end
    nE=[nE(:,1,:) nE(:,2,:)+nE(:,4,:) nE(:,3,:)+nE(:,7,:) nE(:,5,:) nE(:,6,:)+nE(:,8,:) nE(:,9,:)];
    nE=repmat(nE,[1 1 2]);
    nE_tot(:, :, 1)=nE(:, :, 1);
    nE_tot(:, :, 2)=nE(:, :, 2);
    nE_tot(:, :, 3)=nE(:, :, 3);
    nE_tot(:, :, 4)=nE(:, :, 3);
    nE_tot(:, :, 5)=nE(:, :, 1);
    nE_tot(:, :, 6)=nE(:, :, 2); 
    nE=nE_tot;
    al=2*pi*1i*freq.*reshape(L(:, [1 2 3 1 2 3]),[Nt 1 1 6]);
    B=zeros(Nt,1,4,6);

    B(:,1,1,:)=al.*(2+1/3*al.^2);
    
    B(:,1,2:4,:)=repmat(al.^2,[1 1 3 1]).*reshape((p(:,:,[2 3 1 5 6 4])-mean(p,3))./...
        reshape(L(:,[1 2 3 1 2 3]),[Nt 1 6]),[Nt 1 3 6])...
        +repmat(al.^3,[1 1 3 1]).*reshape((p(:,:,[3 1 2 4 5 6])-mean(p,3))./...
        reshape(L(:,[1 2 3 1 2 3]),[Nt 1 6]),[Nt 1 3 6])...
        -repmat(al.^3,[1 1 3 1]).*reshape((p(:,:,[1 2 3 6 4 5])-mean(p,3))./...
        reshape(L(:,[1 2 3 1 2 3]),[Nt 1 6]),[Nt 1 3 6]);

%     B(:,1,1,:)=al.*(2-2*al.^2+4/3*al.^2);
%     
%     B(:,1,2:4,:)=repmat(al.^2,[1 1 3 1]).*reshape((p(:,:,[2 3 1 6 4 5])-mean(p,3))./...
%         reshape(L(:,[1 2 3 1 2 3]),[Nt 1 6]),[Nt 1 3 6])...
%         -repmat(al.^3,[1 1 3 1]).*reshape((p(:,:,[2 3 1 6 4 5])-mean(p,3))./...
%         reshape(L(:,[1 2 3 1 2 3]),[Nt 1 6]),[Nt 1 3 6])...
%         +repmat(al.^3,[1 1 3 1]).*reshape((p(:,:,[3 1 2 5 6 4])-mean(p,3))./...
%         reshape(L(:,[1 2 3 1 2 3]),[Nt 1 6]),[Nt 1 3 6])...
%         -repmat(al.^3,[1 1 3 1]).*reshape((p(:,:,[1 2 3 4 5 6])-mean(p,3))./...
%         reshape(L(:,[1 2 3 1 2 3]),[Nt 1 6]),[Nt 1 3 6]);
    
    A=.5*reshape(nE,[Nt,6,1,6]).*B;
    A(:,[1 4],:,:)=A(:,[1 4],:,:)-A(:,[6],:,:);
    A=A(:,1:5,:,:);
    
    tmp=A;
    %% simple test
%     A=tmp(:,:,:,[1 2 3]);
%     A=tmp;
    %% original TDI
%     A=tmp(:,:,:,[1 2 3]).*reshape(exp(-4*pi*1i*freq*L(:,[3 1 2])) - ones(size(L(:,[1 2 3]))), [Nt, 1, 1, 3])...
%         -tmp(:,:,:,[4 5 6]).*reshape(exp(-4*pi*1i*freq*L(:,[1 2 3])) - ones(size(L(:,[1 2 3]))), [Nt, 1, 1, 3]);
%     A=tmp(:,:,:,[4 5 6]);
%     A=tmp(:,:,:,[1 2 3]).*reshape(exp(-4*pi*1i*freq*L(:,[3 1 2])) - ones(size(L(:,[1 2 3]))), [Nt, 1, 1, 3]);
%     A=tmp(:,:,:,[4 5 6]).*reshape(exp(-4*pi*1i*freq*L(:,[1 2 3])) - ones(size(L(:,[1 2 3]))), [Nt, 1, 1, 3]);

%     A=tmp(:,:,:,[1 2 3]).*reshape(exp(-4*pi*1i*freq*L(:,[3 1 2])) - ones(size(L(:,[1 2 3]))), [Nt, 1, 1, 3])...
%     -tmp(:,:,:,[4 5 6]).*reshape(exp(-4*pi*1i*freq*L(:,[1 2 3])) - ones(size(L(:,[1 2 3]))), [Nt, 1, 1, 3]);
%     A=A./reshape(factor1(:,[1 2 3]),[],1,1,3);
%     AFFT=fftshift(fft(A,[],1),1);
%     AFFT=AFFT.*reshape(factor2(:,[1 2 3]),[],1,1,3);
%     A=ifft(ifftshift(AFFT,1),[],1);
    %% FastGB
    A=tmp(:,:,:,[1 2 3]).*reshape(exp(-4*pi*1i*freq*L(:,[3 1 2])) - ones(size(L(:,[1 2 3]))), [Nt, 1, 1, 3])...
    -tmp(:,:,:,[4 5 6]).*reshape(exp(-4*pi*1i*freq*L(:,[1 2 3])) - ones(size(L(:,[1 2 3]))), [Nt, 1, 1, 3]);
end

