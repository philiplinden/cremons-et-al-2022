function SSA = Hapke(R,p,WLS)
mu=cosd(0);
mu0=cosd(0);
g=0;
h=(-3/8)*log(1-0.4);
B=1/(1+(1/h)*tand(g/2));
for m=1:numel(WLS)
Refl=R(m);
y=@(SSA,x) (SSA./(4*(mu0+mu))).*(p.*(1+B)+(((1-SSA.*mu0.*(((1-sqrt(1-SSA))./(1+sqrt(1-SSA)))+((1-2.*((1-sqrt(1-SSA))./(1+sqrt(1-SSA))).*mu0)/2)*log((1+mu0)./mu0))).^-1).*((1-SSA.*mu.*(((1-sqrt(1-SSA))./(1+sqrt(1-SSA)))+((1-2.*((1-sqrt(1-SSA))./(1+sqrt(1-SSA))).*mu)/2)*log((1+mu)./mu))).^-1))-1);
x=WLS(m);
yx=Refl;
w0=0.5; %Initial Guesses here
OLS=@(SSA) sum((y(SSA,x)-yx).^2); %Ordinary Least Squares Function
opts=optimset('MaxIter',10000,'MaxFunEvals',10000,'TolFun',1E-7,'TolX',1E-7);
SSA(m,1) = fminsearch(OLS, w0, opts);      % Use ‘fminsearch’ to minimise the ‘OLS’ function
end
