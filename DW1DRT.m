%% t = DW1DRT(zs, delta, vel, dep)
% This is a ray tracing program for direct wave in 1-Dimension horizontal layered medium.
% Written by Tche.L. from USTC, 2016.8.
%
% t: a single variable, the calculated travel time for direct wave of the $rlyr$ -th layer with epicenter distance = $delta$.
%
% zs: a single variable, the depth of source point.
% delta: a single variable, the epicenter distance.
% vel: an array whose size is [1, $nlyr$], the velocity of every layer, the kth element is the wave velocity of the kth layer from top to bottom; Endnote: $nlyr$ is the layer number of medium.
% dep: an array whose size is [1, $nlyr$], the depth of every layer, the kth element is the depth of the interface from top to bottom while the 1st one which represents earth surface is zero.
%
% Ref.: 田玥, 陈晓非. 2005. 水平层状介质中的快速两点间射线追踪方法. 地震学报: 27(2), 147-154.
% Example: t = DW1DRT(3500, 5000, vel ,dep) where:
%   vel = [2000, 3200, 1800, 2800, 4000];
%   dep = [0, 1000, 2000, 3050, 3650, 5000];

%%
function t = DW1DRT(zs, delta, vel, dep)
%% Input the parameter.

eps = 1e-6;
itsmax = 10000;

%% Calculate the parameters.

h = diff(dep);
nlyr = length(vel);

if(zs > dep(end))
    error('MyArgumentInvalidError:zs', ...
        'The #1 argument $zs$ is invalid because the source should be above the last layer.');
end
slyr = max(1,find(zs <= dep,1,'first') - 1);

hev(1:slyr - 1) = h(1:slyr - 1);
hev(slyr) = zs - dep(slyr);
% hev(slyr + 1:nlyr) = 0;

[vmax, imax] = max(vel(1:slyr));
hmax = h(imax);
ek = vel(1:slyr)./vmax;

a = sum(ek.*hev)/hmax;
b = sum(ek(1:imax - 1).*hev(1:imax - 1)./sqrt(1 - ek(1:imax - 1).^2)) ...
    + sum(ek(imax + 1:slyr).*hev(imax + 1:slyr)./sqrt(1 - ek(imax + 1:slyr).^2));
DltC = a*b/(a - 1);

%% Calculate the p value.

if(zs == 0)
    p = 1/vel(1);
else
    if(delta < DltC)
        q = delta/a;
    else
        q = delta - b;
    end
    its = 0;
    fq = q*sum(ek.*hev./sqrt(hmax^2 + (1 - ek.^2).*q^2)) - delta;
    while(abs(fq/delta) > eps && its <= itsmax)
        its = its + 1;
        q = q - fq/(hmax^2*sum(ek.*hev./(sqrt(hmax^2 + (1 - ek.^2).*q^2)).^3));
        fq = q*sum(ek.*hev./sqrt(hmax^2 + (1 - ek.^2).*q^2)) - delta;
    end
    if(its > itsmax)
        warning('WARNING:MyMaxIteratingTimes',['The iterating time of Newton method is beyond the maximum time of iterations and', ...
            ' the iterating procedure may be divergent, so the result may be not right!']);
    end
    p = q/(vmax*sqrt(hmax^2 + q^2));
end

%% Calculate the travel time.

eta = sqrt(1./(vel(1:slyr).^2) - p^2);
t = p*delta + sum(hev.*eta);

%% Calculate the ray path.

dp = NaN*ones(slyr + 1,2);
dp(1,:) = [0, zs];

if(zs == 0)
    dp(2,:) = [delta, 0];
else
    for i = 2:1:(slyr + 1)
        lyrnum = slyr - i + 2;
        dp(i,2) = dep(lyrnum);
        theta = asin(p*vel(lyrnum));
        dp(i,1) = dp(i - 1,1) + hev(lyrnum)*tan(theta);
    end
end

%% Plot the ray path.

figure; hold on;
plot(dp(:,1),dp(:,2),'-b'); axis ij;% axis equal;
plot(dp(1,1),dp(1,2),'rp','MarkerFaceColor','r');
plot(dp(end,1),dp(end,2),'k^','MarkerFaceColor','k');
for i = 1:1:(nlyr + 1)
    plot([0 delta],[dep(i) dep(i)],'--k');
end
xlabel('The epicenter distance'); ylabel('The depth'); hold off;
title({sprintf('The ray path with epicenter distance = %.2f for direct wave',delta), ...
    sprintf('And the travel time is: %.2f s',t)});

end


