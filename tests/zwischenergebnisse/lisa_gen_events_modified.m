function [Y, dApsi, phaseShift, xi, ep, ec]=lisa_gen_events_modified(t,pr,ps,n,L,parameters,freq0)
pr=pr./299792458;
ps=ps./299792458;
L=L./299792458;

M=size(parameters,1);
Yall=zeros(length(t),6,M);
phaseShiftAll=zeros(length(t),6,M);
dApsiAll=zeros(length(t),6,M);
xiAll=zeros(length(t),6,M);

for m=1:M
    lambda=parameters(m,1);
    beta=parameters(m,2);
    psi=parameters(m,3);
    amp=parameters(m,4);
    iota=parameters(m,5);
    freq=parameters(m,6);
    phi0=-1*parameters(m,7);
    freqdot=parameters(m,8);
    
    k=zeros(3,1);
    k(1)=-cos(beta)*cos(lambda);
    k(2)=-cos(beta)*sin(lambda);
    k(3)=-sin(beta);
    
    [Ep, Ec, ep ,ec]=get_polarization_tensors(lambda, beta, psi);
    
    y=zeros(length(t),6);
    dApsi=zeros(length(t),6);
    phaseShift=zeros(length(t),6);
    yRecon=zeros(length(t),6);
    xi=zeros(length(t),6);
    
   for jord=1:6
        [hp_r,hc_r]=h_components(t,-pr(:,:,jord)*k(:),iota,amp,freq,phi0,freqdot,freq0(m));
        [hp_s,hc_s]=h_components(t,-L(:,jord)-ps(:,:,jord)*k(:),iota,amp,freq,phi0,freqdot,freq0(m));
        
        psir=.5*(sum((n(:,:,jord)*Ep).*n(:,:,jord),2).*hp_r+sum((n(:,:,jord)*Ec).*n(:,:,jord),2).*hc_r)./(1-n(:,:,jord)*k(:));
        psis=.5*(sum((n(:,:,jord)*Ep).*n(:,:,jord),2).*hp_s+sum((n(:,:,jord)*Ec).*n(:,:,jord),2).*hc_s)./(1-n(:,:,jord)*k(:));

        yRecon(:,jord)=-1*(psis-psir); 
        yRecon(:,jord)=yRecon(:,jord)./(4*pi*1i*(freq+freqdot.*(t-ps(:,:,jord)*k(:)-L(:,jord))).*L(:,jord));
        
        dApsi(:,jord)=(sum((n(:,:,jord)*Ep).*n(:,:,jord),2).*amp*(1+cos(iota)*cos(iota))+...
            (2*1i*amp*cos(iota))*sum((n(:,:,jord)*Ec).*n(:,:,jord),2));

        omega=freq+freqdot.*(t-ps(:,:,jord)*k(:));
        phaseShift(:, jord)=phaseShift_fun(t, -ps(:,:,jord)*k(:), freq, freqdot, phi0, freq0(m));
        phaseShift(:, jord)=phaseShift(:, jord).*(exp(2*pi*1i*omega.*L(:,jord).*(1+n(:,:,jord)*k(:)))-1);
        phaseShift(:,jord)=0.5.*phaseShift(:,jord)./(1+n(:,:,jord)*k(:));
        phaseShift(:, jord)=phaseShift(:, jord)./(4*pi*1i*(freq+freqdot.*(t-ps(:,:,jord)*k(:))).*L(:,jord));
        
%         omega=freq+freqdot.*(t-ps(:,:,jord)*k(:)-L(:,jord));
%         phaseShift(:, jord)=phaseShift_fun(t, -L(:,jord)-ps(:,:,jord)*k(:), freq, freqdot, phi0, freq0(m));
%         phaseShift(:, jord)=phaseShift(:, jord).*(exp(2*pi*1i*omega.*L(:,jord).*(1-n(:,:,jord)*k(:)))-1);
%         phaseShift(:,jord)=0.5.*phaseShift(:,jord)./(1-n(:,:,jord)*k(:));
%         phaseShift(:, jord)=phaseShift(:, jord)./(4*pi*1i*(freq+freqdot.*(t-ps(:,:,jord)*k(:)-L(:,jord))).*L(:,jord));
        y(:,jord)=dApsi(:,jord).*phaseShift(:,jord);
        
        xi(:, jord)=(t-ps(:,:,jord)*k(:));
        
        if m==1
            err=norm(y(:,jord)-yRecon(:,jord),'fro')/norm(y(:,jord),'fro')

            clf
            subplot(2,2,1)
            plot(real(y(:,jord))); hold on
            plot(real(yRecon(:,jord)),'--')
            title(num2str(err))

            subplot(2,2,2)
            plot(imag(y(:,jord))); hold on
            plot(imag(yRecon(:,jord)),'--')

            subplot(2,2,3)
            yFFT=fftshift(fft(fftshift(y(:,jord))));
            yReconFFT=fftshift(fft(fftshift(yRecon(:,jord))));
            plot(real(yFFT)); hold on
            plot(real(yReconFFT),'--')
            xlim([numel(t)/2-25 numel(t)/2+25])

            subplot(2,2,4)
            plot(imag(yFFT)); hold on
            plot(imag(yReconFFT),'--')
            xlim([numel(t)/2-25 numel(t)/2+25])

            drawnow
        end
   end
    Yall(:,:,m)=y;
    dApsiAll(:,:,m)=dApsi;
    phaseShiftAll(:,:,m)=phaseShift;
    xiAll(:,:,m)=xi;
end
Y=Yall;
phaseShift=phaseShiftAll;
dApsi=dApsiAll;
xi=xiAll;
end


function [Ep, Ec, ep ,ec]=get_polarization_tensors(la, be, psi)
    ep=-1.*[sin(be)^2*cos(la)^2-sin(la)^2 0.5*sin(2*la)*(1+sin(be)^2) -0.5*sin(2*be)*cos(la); ...
        0.5*sin(2*la)*(1+sin(be)^2) sin(be)^2*sin(la)^2-cos(la)^2 -0.5*sin(2*be)*sin(la);...
        -0.5*sin(2*be)*cos(la) -0.5*sin(2*be)*sin(la) cos(be)^2];
    ec=[sin(be)*sin(2*la) -cos(2*la)*sin(be) -cos(be)*sin(la); ...
        -cos(2*la)*sin(be) -sin(be)*sin(2*la) cos(be)*cos(la); ...
        -cos(be)*sin(la) cos(be)*cos(la) 0];
    
    Ep=cos(2*psi).*ep-sin(2*psi).*ec;
    Ec=sin(2*psi).*ep+cos(2*psi).*ec;
end

function ps=phaseShift_fun(t, ts, freq, freqdot, phi0, freq0)
    ps=exp(2*pi*1i*(freq*ts+(freq-freq0)*t)+1i*phi0+pi*1i*freqdot*(ts+t).^2);
end


function [hp,hc]=h_components(t,ts,iota,amp,freq,phi0,freqdot, freq0)
    hp=amp*(1+cos(iota)*cos(iota))*exp(2*pi*1i*(freq*ts+(freq-freq0)*t)+1i*phi0+pi*1i*freqdot*(ts+t).^2); 
    hc=2*amp*cos(iota)*(1i)*exp(2*pi*1i*(freq*ts+(freq-freq0)*t)+1i*phi0+pi*1i*freqdot*(ts+t).^2);
end
