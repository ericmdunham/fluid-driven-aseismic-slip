% script to reproduce plots in
% Dunham, E. M., Fluid-driven aseismic fault slip with permeability enhancement and dilatancy, submitted.

savePlots = false; % set to true to save .eps figures

epsilon = [0 1 10 100]; M = length(epsilon);
lambda = logspace(-2,1,10000)'; N = length(lambda);
T2 = zeros(N,M); T2_BV = zeros(N,1);
T3 = zeros(N,M); T2_SL = zeros(N,1);
T2lambdaSmall = zeros(N,M); T3lambdaSmall = zeros(N,M);

for m=1:M
    for n=1:N
        T2(n,m) = fluidDrivenAseismicSlip.evalT2(lambda(n),epsilon(m));
        T3(n,m) = fluidDrivenAseismicSlip.evalT3(lambda(n),epsilon(m));
        if m==1, T2_BV(n,1) = fluidDrivenAseismicSlip.evalT2_BV(lambda(n)); end
        if m==1, T3_SL(n,1) = fluidDrivenAseismicSlip.evalT3_SL(lambda(n)); end
        %if T(n,m)<0, T(n,m)=0; end
        if isnan(T2(n,m)), T2(n,m)=-n; end % need value for interpolation
        if isnan(T3(n,m)), T3(n,m)=-n; end % need value for interpolation
        if isinf(T2(n,m)), T2(n,m)=-n; end % need value for interpolation
        if isinf(T3(n,m)), T3(n,m)=-n; end % need value for interpolation
        T2lambdaSmall(n,m) = 1-4/pi*(1+epsilon(m))*lambda(n)^2; % asymptotics
        T3lambdaSmall(n,m) = 1/lambda(n)^2+1-epsilon(m)-log(4); % asymptotics
    end
end

figure(1),clf
semilogx(lambda,T2,lambda,T2_BV,'--',lambda,T2lambdaSmall,'k--')
xlabel('\lambda = a(t)/(4\alphat)^{1/2}')
ylabel('stress-injection parameter, T')
ylim([0 1])
lg=legend(strcat('\epsilon = ',num2str(epsilon')));
str = lg.String; str{M+1} = 'BV19'; legend(str)
title('2D constant pressure injection')
if savePlots,print -depsc2 Tlambda2D, end

figure(2),clf
loglog(lambda,T3,lambda,T3_SL,'--',lambda,T3lambdaSmall,'k--')
xlabel('\lambda = a(t)/(4\alphat)^{1/2}')
ylabel('stress-injection parameter, T')
ylim([1e-1 1e2])
lg=legend(strcat('\epsilon = ',num2str(epsilon')));
str = lg.String; str{M+1} = 'SLBV22'; legend(str)
title('3D constant rate injection')
if savePlots,print -depsc2 Tlambda3D, end

% nondimensional pressure
xi = [linspace(0,0.01,10) linspace(0.01,1,100) linspace(1,2,100)]'; I = 1:110;
Tplot2D = 0.5;
for m=1:M
    la2(m) = interp1(T2(:,m),lambda,Tplot2D);
    g2(:,m) = fluidDrivenAseismicSlip.evalG2(xi,la2(m),epsilon(m));
end
g2(I(end)+1,:) = 0;
la2_BV = interp1(T2_BV,lambda,Tplot2D);
g2_BV = fluidDrivenAseismicSlip.evalG2_BV(xi,la2_BV);
g2_BV(I(end)+1) = 0;

Tplot3D = 0.5;
for m=1:M
    la3(m) = interp1(T3(:,m),lambda,Tplot3D);
    g3(:,m) = fluidDrivenAseismicSlip.evalG3(xi,la3(m),epsilon(m));
end
g3(I(end)+1,:) = 0;
la3_SL = interp1(T3_SL,lambda,Tplot3D);
g3_SL = fluidDrivenAseismicSlip.evalG3_SL(xi,la3_SL);
g3_SL(I(end)+1) = 0;

figure(3)
plot(xi,g2,xi,g2_BV,'--')
xlabel('\xi = x / \lambda(4\alphat)^{1/2}')
ylabel('pressure change, p/\Deltap = g(\xi)')
legend(num2str(epsilon'))
lg=legend(strcat('\epsilon = ',num2str(epsilon')));
str = lg.String; str{M+1} = 'BV19'; legend(str)
title(['2D constant pressure injection, T = ' num2str(Tplot2D)])
if savePlots,print -depsc2 p2D, end

figure(4)
plot(xi,g3,xi,g3_SL,'--')
xlabel('\xi = r / \lambda(4\alphat)^{1/2}')
ylabel('pressure change, p/\Deltap = g(\xi)')
legend(num2str(epsilon'))
lg=legend(strcat('\epsilon = ',num2str(epsilon')));
str = lg.String; str{M+1} = 'SLBV22'; legend(str)
title(['3D constant rate injection, T = ' num2str(Tplot3D)])
ylim([-0.5 2])
if savePlots,print -depsc2 p3D, end

% nondimensional stress drop

figure(5)
plot(xi(I),g2(I,:)-Tplot2D,xi(I),g2_BV(I)-Tplot2D,'--')
line([0 1],[0 0],'color','k','linestyle','--')
xlabel('\xi = x / \lambda(4\alphat)^{1/2}')
ylabel('stress drop, \Delta\tau/(f\Deltap) = g(\xi)-T')
legend(num2str(epsilon'))
lg=legend(strcat('\epsilon = ',num2str(epsilon')));
str = lg.String; str{M+1} = 'BV19'; legend(str)
title(['2D constant pressure injection, T = ' num2str(Tplot2D)])
if savePlots,print -depsc2 tau2D, end

figure(6)
plot(xi(I),g3(I,:)-Tplot3D,xi(I),g3_SL(I)-Tplot3D,'--')
line([0 1],[0 0],'color','k','linestyle','--')
xlabel('\xi = r / \lambda(4\alphat)^{1/2}')
ylabel('stress drop, \Delta\tau/(f\Deltap) = g(\xi)-T')
legend(num2str(epsilon'))
lg=legend(strcat('\epsilon = ',num2str(epsilon')));
str = lg.String; str{M+1} = 'SLBV22'; legend(str)
title(['3D constant rate injection, T = ' num2str(Tplot3D)])
if savePlots,print -depsc2 tau3D, end

% repeat 3D for larger T

Tplot3D = 30;
for m=1:M
    la3(m) = interp1(T3(:,m),lambda,Tplot3D);
    g3(:,m) = fluidDrivenAseismicSlip.evalG3(xi,la3(m),epsilon(m));
end
g3(I(end)+1,:) = 0;

figure(7)
plot(xi,g3)
xlabel('\xi = r / \lambda(4\alphat)^{1/2}')
ylabel('pressure change, p/\Deltap = g(\xi)')
legend(num2str(epsilon'))
lg=legend(strcat('\epsilon = ',num2str(epsilon')));
title(['3D constant rate injection, T = ' num2str(Tplot3D)])
xlim([0 1.2])
if savePlots,print -depsc2 p3D_largeT, end

% nondimensional stress drop

figure(8)
plot(xi(I),g3(I,:)-Tplot3D)
line([0 1],[0 0],'color','k','linestyle','--')
xlabel('\xi = r / \lambda(4\alphat)^{1/2}')
ylabel('stress drop, \Delta\tau/(f\Deltap) = g(\xi)-T')
legend(num2str(epsilon'))
lg=legend(strcat('\epsilon = ',num2str(epsilon')));
title(['3D constant rate injection, T = ' num2str(Tplot3D)])
if savePlots,print -depsc2 tau3D_largeT, end

% repeat 3D for different T, all epsilon=0

T = [1 3 10 30]; nT = length(T);
for n=1:nT
    la3T(n) = interp1(T3(:,1),lambda,T(n));
    g3T(:,n) = fluidDrivenAseismicSlip.evalG3(xi,la3T(n),0);
    g3TsmallLambda(:,n) = fluidDrivenAseismicSlip.evalG3smallLambda(xi,la3T(n),0);
end
g3T(I(end)+1,:) = 0;
g3TsmallLambda(I(end)+1,:) = 0;

figure(9)
for n=1:nT
    plot(xi,g3T(:,n)/T(n))
    if n==1,hold on,end
end
for n=1:nT
    plot(xi,g3TsmallLambda(:,n)/T(n),'k--')
end
hold off
xlabel('\xi = r / \lambda(4\alphat)^{1/2}')
ylabel('pressure change, p/(T\Deltap) = g(\xi)/T')
legend(num2str(T'))
lg=legend(strcat('T = ',num2str(T')));
title('3D constant rate injection, \epsilon = 0')
ylim([0 5]),xlim([0 1.2])
if savePlots,print -depsc2 p3D_varyT, end

figure(10)
for n=1:nT
    plot([xi*la3T(n); 1.2],[g3T(:,n); 0])
    if n==1,hold on,end
end
plot(xi*la3_SL,g3_SL,'--')
for n=1:nT
    plot([xi*la3T(n); 1.2],[g3TsmallLambda(:,n); 0],'k--')
end
hold off
xlabel('r/(4\alphat)^{1/2}')
ylabel('pressure change, p/\Deltap = g(\xi)')
legend(num2str(T'))
lg=legend(strcat('T = ',num2str(T')));
str = lg.String; str{nT+1} = 'SLBV22'; legend(str)
title('3D constant rate injection, \epsilon = 0')
ylim([0 40]),xlim([0 1.2])
if savePlots,print -depsc2 p3D_varyT3_alt, end

