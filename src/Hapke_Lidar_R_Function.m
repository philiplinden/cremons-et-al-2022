function R = Hapke(w,p,WLS)
mu=cos(0);
mu0=cos(0);
g=0;
SSATot=w;
for m=1:numel(SSATot)
    if SSATot(m)>1
        SSATot(m)=1;
    else
    end
end
h=(-3/8)*log(1-0.41); %for lunar regolith from Li et al 2011
B=1/(1+(1/h)*tand(g/2));
gamma=sqrt(1-SSATot);
r0=(1-gamma)./(1+gamma);
H0=(1-SSATot.*mu0.*(r0+((1-2.*r0.*mu0)/2)*log((1+mu0)./mu0))).^-1;
H=(1-SSATot.*mu.*(r0+((1-2.*r0.*mu)/2)*log((1+mu)./mu))).^-1;
R=(SSATot./(4*(mu0+mu))).*(p.*(1+B)+(H0.*H)-1); %Radiance Coefficeint