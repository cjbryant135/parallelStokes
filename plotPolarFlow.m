clear all;
close all;

set(0,'defaulttextinterpreter','latex');

outDir = './';

suffix = '.out';

uPOut = fopen([outDir 'u' suffix]);
vPOut = fopen([outDir 'v' suffix]);
pPOut = fopen([outDir 'p' suffix]);

%read polar data
upmaxi = fread(uPOut, 1, 'int');
upmaxj = fread(uPOut, 1, 'int');
up     = fread(uPOut, [upmaxj, upmaxi], 'double');
upx    = fread(uPOut, [upmaxj, upmaxi], 'double');
upy    = fread(uPOut, [upmaxj, upmaxi], 'double');

vpmaxi = fread(vPOut, 1, 'int');
vpmaxj = fread(vPOut, 1, 'int');
vp     = fread(vPOut, [vpmaxj, vpmaxi], 'double');
vpx    = fread(vPOut, [vpmaxj, vpmaxi], 'double');
vpy    = fread(vPOut, [vpmaxj, vpmaxi], 'double');

ppmaxi = fread(pPOut, 1, 'int');
ppmaxj = fread(pPOut, 1, 'int');
pp     = fread(pPOut, [ppmaxj, ppmaxi], 'double');
ppx    = fread(pPOut, [ppmaxj, ppmaxi], 'double');
ppy    = fread(pPOut, [ppmaxj, ppmaxi], 'double');

figure(1)
subplot(1,3,1)
contourf(upx,upy,up,100)
colorbar
axis square
set(gca,'fontsize',18)
title('$u$')

subplot(1,3,2)
contourf(vpx,vpy,vp,100)
colorbar
axis square
set(gca,'fontsize',18)
title('$v$')

subplot(1,3,3)
contourf(ppx,ppy,pp,100)
colorbar
axis square
set(gca,'fontsize',18)
title('$p$')


figure(2)
subplot(1,3,1)
contourf(upx,upy,up-upx.*sqrt(upx.^2+upy.^2),100)
colorbar
axis square
set(gca,'fontsize',18)
title('$u$ Errors')

subplot(1,3,2)
contourf(vpx,vpy,vp--3*vpy.*sqrt(vpx.^2+vpy.^2),100)
colorbar
axis square
set(gca,'fontsize',18)
title(['$v$ Errors'])

subplot(1,3,3)
contourf(ppx,ppy,pp-ppy,100)
colorbar
axis square
set(gca,'fontsize',18)
title('$p$ Errors')



rIJ = sqrt(ppx.^2+ppy.^2);
tIJ = atan2(ppy,ppx);

dr = rIJ(1,2)-rIJ(1,1);
dth = tIJ(2,1)-tIJ(1,1);

%pad theta-velocity
vp(end+1,:) = vp(1,:);
pDiv = poleDiv(up,vp,dr,dth,rIJ(1,:));

figure(3)
contourf(ppx,ppy,pDiv,100)
colorbar
axis square
set(gca,'fontsize',18)
title('$\nabla \cdot u$')

function d = poleDiv(u,v,dr,dth,r)
%compute divergence of polar flow - r should be a vector of rC values
    d = (u(:,2:end)-u(:,1:end-1))./dr + (u(:,2:end)+u(:,1:end-1))./r./2 + ...
        (v(2:end,:)-v(1:end-1,:))./r./dth;
end



