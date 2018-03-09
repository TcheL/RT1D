%% t = RW1DRT(zs, delta, rlyr, vel, dep)
% This is a ray tracing program for reflected wave in 1-Dimension horizontal layered medium.
% Written by Tche.L. from USTC, 2016.4.
%
% t: a single variable, the calculated travel time for reflected wave of the $rlyr$ -th layer with epicenter distance = $delta$ .
%
% zs: a single variable, the depth of source point.
% delta: a single variable, the epicenter distance.
% rlyr: a single variable, the layer number of reflected layer.
% vel: a vector whose size is [1, $nlyr$], the velocity of every layer, the kth element is the wave velocity of the kth layer from top to bottom; Endnote: $nlyr$ is the layer number of medium.
% dep: a vector whose size is [1, $nlyr$], the depth of every layer, the kth element is the depth of the kth interface from top to bottom while the 1st which represents earth surface is zero).
%
% Ref.: 田玥, 陈晓非. 2005. 水平层状介质中的快速两点间射线追踪方法. 地震学报: 27(2), 147-154.
% Example: t = RW1DRT(0, 5000, 5, vel ,dep) where:
%   vel = [2000, 3200, 1800, 2800, 4000];
%   dep = [0, 1000, 2000, 3050, 3650, 5000];

%%
function t = RW1DRT(zs, delta, rlyr, vel, dep)
%% Input the parameter.

eps = 1e-6;                                                                                         % precision
itsmax = 10000;                                                                                     % the maximum number of iterations of Newton method.

%% Calculate the parameters.

h = diff(dep);                                                                                      % thickness
nlyr = length(vel);

if(zs > dep(end))
    error('MyArgumentInvalidError:zs', ...
        'The #1 argument $zs$ is invalid because the source should be above the last layer.');
end
slyr = max(1,find(zs <= dep,1,'first') - 1);
if(rlyr < slyr)
    error('MyArgumentInvalidError:rlyr', ...
        'The #3 argument $rlyr$ is invalid because the reflected layer must be below source point.');
end

hev(1:slyr - 1) = h(1:slyr - 1);
hev(slyr) = dep(slyr + 1) - zs + h(slyr);
hev(slyr + 1:rlyr) = 2.*h(slyr + 1:rlyr);

[vmax, imax] = max(vel(1:rlyr));
hmax = h(imax);
ek = vel(1:rlyr)./vmax;

a = sum(ek.*hev)/hmax;
b = sum(ek(1:imax - 1).*hev(1:imax - 1)./sqrt(1 - ek(1:imax - 1).^2)) ...
    + sum(ek(imax + 1:rlyr).*hev(imax + 1:rlyr)./sqrt(1 - ek(imax + 1:rlyr).^2));
DltC = a*b/(a - 1);

%% Calculate the p value.

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

%% Calculate the travel time.

eta = sqrt(1./(vel(1:rlyr).^2) - p^2);
t = p*delta + sum(hev.*eta);

%% Calculate the ray path.

rp = NaN*ones(2*rlyr - slyr + 2,2);
rp(1,:) = [0, zs];
rp(2,2) = dep(slyr + 1);
theta = asin(p*vel(slyr));
rp(2,1) = (rp(2,2) - zs)*tan(theta);

for i = 3:1:(rlyr - slyr + 2)
    lyrnum = slyr + i - 2;
    rp(i,2) = rp(i - 1,2) + h(lyrnum);
    theta = asin(p*vel(lyrnum));
    rp(i,1) = rp(i - 1,1) + h(lyrnum)*tan(theta);
end
for i = (rlyr - slyr + 3):1:(2*rlyr - slyr + 2)
    lyrnum = 2*rlyr - slyr - i + 3;
    rp(i,2) = rp(i - 1,2) - h(lyrnum);
    theta = asin(p*vel(lyrnum));
    rp(i,1) = rp(i - 1,1) + h(lyrnum)*tan(theta);
end

%% Plot the ray path.

figure; hold on;
plot(rp(:,1),rp(:,2),'-b'); axis ij;% axis equal;
plot(rp(1,1),rp(1,2),'rp','MarkerFaceColor','r');
plot(rp(end,1),rp(end,2),'k^','MarkerFaceColor','k');
for i = 1:1:(nlyr + 1)
    plot([0 delta],[dep(i) dep(i)],'--k');
end
hold off;
xlabel('The epicenter distance'); ylabel('The depth');
title({sprintf('The ray path with epicenter distance = %.2f for reflected wave of the %dth layer',delta,rlyr), ...
    sprintf('And the travel time is: %.2f s',t)});

end