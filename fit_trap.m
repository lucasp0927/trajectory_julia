function [ prob, xstart, xend, ystart, yend] = fit_trap( input, x_temp, y_temp )
xin = squeeze(input(1,:,:));
yin = squeeze(input(2,:,:));
zin = squeeze(input(3,:,:));
[xout,yout,zout]=prepareSurfaceData(xin,yin,zin);
[fitr,gof]=createFit(xout,yout,zout);
xx = squeeze(input(1,:,:));
yy = squeeze(input(2,:,:));
fitout = fitr(xx,yy);
trap = squeeze(input(3,:,:));
residual = trap-fitout;
% figure;
% subplot(3,1,1);
% h=pcolor(yy,xx,trap);
% set(h,'edgecolor','none');
% title('trap potential');
% xlabel('nm');
% ylabel('nm');
% colorbar();
% subplot(3,1,2);
% h=pcolor(yy,xx,fitout);
% set(h,'edgecolor','none');
% title('fit potential');
% xlabel('nm');
% ylabel('nm');
% colorbar();
% subplot(3,1,3);
% h=pcolor(yy,xx,residual);
% title('residual');
% set(h,'edgecolor','none');
% xlabel('nm');
% ylabel('nm');
% colorbar();

%min_u = min(min(trap));
KB = 1.38065e-17;
x_temperature = x_temp;
y_temperature = y_temp;
%coeffs = coeffvalues(fitr);
%disp(coeffs)
%todo
xmean = mean(xout);
xstd = std(xout);
ymean = mean(yout);
ystd = std(yout);
xx_old = xx;
yy_old = yy;
xx = (xx-xmean)/xstd;
yy = (yy-ymean)/ystd;

FX = @(x) fitr.p10*x+fitr.p20*(x.^2)+fitr.p30*(x.^3)+fitr.p40*(x.^4);
FY = @(y) fitr.p01*y+fitr.p02*(y.^2)+fitr.p03*(y.^3)+fitr.p04*(y.^4);
FXY = @(x,y) FX(x)+FY(y);
%UX = FX(xx);
%UY = FY(yy);
%figure
%imagesc(UX)
%title('UX')
%figure
%imagesc(UY)
%title('UY')
p = exp(-FX(xx)/x_temperature).*exp(-FY(yy)/y_temperature).*exp(-residual/y_temperature);
prob = p/((max(max(p)))*1.01);

xstart = xx_old(1,1);
xend = xx_old(end,1);
ystart = yy_old(1,1);
yend = yy_old(1,end);
%figure;
%h=pcolor(yy_old,xx_old,p);
%set(h,'edgecolor','none');
%coeff = [fitr.p00, fitr.p10, fitr.p01, fitr.p20, fitr.p11, fitr.p02, fitr.p30, fitr.p21, fitr.p12, fitr.p03, fitr.p40, fitr.p31, fitr.p22, fitr.p13, fitr.p04];
end
