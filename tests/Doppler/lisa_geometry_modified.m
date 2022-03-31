function [pr,ps,n,L]=lisa_geometry_modified(t)
    Nt=numel(t);
%     L0=2500000000/299792458;
    L0=2500000000;
    ps=zeros(Nt,3,6);pr=zeros(Nt,3,6);
    pord=[1 2;2 3;3 1;2 1;3 2;1 3];Nw = size(pord);Nw=Nw(1);
    L=zeros(Nt,Nw);n=zeros(Nt,3,Nw);
    for jord=1:Nw
      L(:,jord)=L0;
      pr(:,:,jord)=get_p(t,pord(jord,2));
      ps(:,:,jord)=get_p(t,pord(jord,1));
      n(:,:,jord)=pr(:,:,jord)-ps(:,:,jord);
      n(:,:,jord)=n(:,:,jord)/L0;
    end
end

function p=get_p(t,s)
    XI_0=3*pi/2;
    ETA_0 = 0.0;
%     OMEGA = (2.0*pi)/(365.256363004*24*60*60);
    OMEGA=2.0*pi*3.168753575e-8;
    L = 2500000000;
%     R=149597870700;
    R=1.49597870660e11;
%     L = 2500000000/299792458;
%     R=149597870700/299792458;
%     ECC = L/(2*sqrt(3)*R); 
    ECC=0.0048241852;
    beta(1)=ETA_0+XI_0-3*pi/2;
    beta(2)=ETA_0+XI_0-(3*pi/2-2*pi/3);
    beta(3)=ETA_0+XI_0-(3*pi/2-4*pi/3);
    alpha=OMEGA*t+ETA_0;
    p(:, 1) = R*(cos(alpha)+ECC*(sin(alpha).*cos(alpha)*sin(beta(s)) - (1+sin(alpha).*sin(alpha))*cos(beta(s))));
    p(:, 2) = R*(sin(alpha)+ECC*(sin(alpha).*cos(alpha)*cos(beta(s)) - (1+cos(alpha).*cos(alpha))*sin(beta(s))));
    p(:, 3) = -R*sqrt(3.0)*ECC*cos(alpha-beta(s));
end
