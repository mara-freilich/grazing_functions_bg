%% parameters and options
cols = brewermap(3,'Dark2');
load('Timeseries_AFTER_processing.mat')
xdata = day(TimeseriesAFTERprocessing.date,'doy'); xdata(xdata < xdata(1)) = xdata(xdata < xdata(1))+365;
ydata = smoothdata(TimeseriesAFTERprocessing.accumulationratemedian,'sgolay');

rmvstart = 14; rmvend = 29; xdata1 = xdata; ydata1 = ydata;
ydata1(rmvstart:rmvend) = interp1([xdata(rmvstart) xdata(rmvend)],[ydata(rmvstart) ydata(rmvend)],xdata(rmvstart:rmvend));

[pmlII,zmlII,pmlIII,zmlIII,dpdtII,dpdtIII,BII,BIII] = run_models;

t = (230:595)';
[mld,~,~,~,~,~] = mldmodel(t);

% compute the equilibrium concentration over one annual cyclee
[PIII,ZIII] = eq_typeIII(BIII);
[PII,ZII] = eq_typeII(BII);

%% compare observations and models (Figure 3)
figure
subplot(3,2,[1 3 5])
yIII = dpdtIII(t+365*8); yII = dpdtII(t+365*8);
plot(t,0*t,'k'); hold on
plot(t,yIII,'color',cols(1,:),'linewidth',2)
plot(t,yII,'color',cols(2,:),'linewidth',2)
plot(xdata,ydata,'k','linewidth',2)
plot(xdata,ydata1,'k','linewidth',1)
plot(xdata,TimeseriesAFTERprocessing.accumulationratefirstquartile,'color',[0.3 0.3 0.3])
plot(xdata,TimeseriesAFTERprocessing.accumulationratethirdquartile,'color',[0.3 0.3 0.3])
grid on
hold off
x1 = 230:30:600;
set(gca,'fontsize',15,'xtick',x1,'xticklabels',mod(x1,365))
xlim([230 595])
xlabel('day of year')
ylabel('net population growth rate (d^{-1})')

subplot(3,2,[2 4])
plot(t,pmlIII(t+365*8),'color',cols(1,:),'linewidth',2)
hold on
plot(t,pmlII(t+365*8),'color',cols(2,:),'linewidth',2)

plot(xdata,TimeseriesAFTERprocessing.surfacephytoplanktonmedian,'k','linewidth',2);
hold on
plot(xdata,TimeseriesAFTERprocessing.surfacephytoplanktonfirstquartile,'color',[0.3 0.3 0.3])
plot(xdata,TimeseriesAFTERprocessing.surfacephytoplanktonthirdquartile,'color',[0.3 0.3 0.3])
grid on
hold off
x1 = 230:30:600;
set(gca,'fontsize',15,'xtick',x1,'xticklabels','')
xlim([230 595])
ylabel({'surface phytoplankton';'concentration (mg C m^{-3})'})
legend('type III','type II','observations')

subplot(3,2,6)
plot(t,-mld,'color',cols(1,:),'linewidth',2)
hold on
plot(t,-mld,'color',cols(2,:),'linestyle','--','linewidth',2)
plot(xdata,TimeseriesAFTERprocessing.mldmedian,'k','linewidth',2)
plot(xdata,TimeseriesAFTERprocessing.mldfirstquartile,'color',[0.3 0.3 0.3])
plot(xdata,TimeseriesAFTERprocessing.mldthirdquartile,'color',[0.3 0.3 0.3])
grid on
hold off
x1 = 230:30:600;
set(gca,'fontsize',15,'xtick',x1,'xticklabels',mod(x1,365))
xlim([230 595])
ylabel('mixed layer depth (m)')
xlabel('day of year')

fig = figure;
plot(t,BIII(2)*zmlIII(t+365*8).*(pmlIII(t+365*8)./(pmlIII(t+365*8).^2+BIII(3)^2)),'color',cols(1,:),'linewidth',2)
ylabel('grazing rate (type III) (day^{-1})')
hold on
plot(t,BII(2)*zmlII(t+365*8)./(pmlII(t+365*8)+BII(3)),'color',cols(2,:),'linewidth',2)
ylabel('grazing rate (day^{-1})')
grid on
x1 = xticks;
set(gca,'fontsize',20,'xtick',x1,'xticklabels',mod(x1,365))
xlim([230 595])
xlabel('day of year')

%% phase portrait (Figure 4)
figure
plot(pmlII(t+365*8),zmlII(t+365*8),'color',cols(2,:),'linewidth',2)
hold on
plot(pmlIII(t+365*8),zmlIII(t+365*8),'color',cols(1,:),'linewidth',2)
plot(PII,ZII,'color',[0.5 0.5 0.5],'linewidth',2)
plot(PIII,ZIII,'color',[0.5 0.5 0.5],'linewidth',2)
scatter(pmlII(365*8),zmlII(365*8),50,cols(2,:),'filled'); text(pmlII(365*8)+0.5,zmlII(365*8),'0')
scatter(pmlIII(365*8),zmlIII(365*8),50,cols(1,:),'filled'); text(pmlIII(365*8)+0.5,zmlIII(365*8),'0')
doy = 91;
scatter(pmlII(doy+365*8),zmlII(doy+365*8),50,cols(2,:),'filled'); text(pmlII(doy+365*8)+0.5,zmlII(doy+365*8),'91')
scatter(pmlIII(doy+365*8),zmlIII(doy+365*8),50,cols(1,:),'filled'); text(pmlIII(doy+365*8)+0.5,zmlIII(doy+365*8),'91')
doy = 138;
scatter(pmlII(doy+365*8),zmlII(doy+365*8),50,cols(2,:),'filled'); text(pmlII(doy+365*8)+0.5,zmlII(doy+365*8),num2str(doy))
scatter(pmlIII(doy+365*8),zmlIII(doy+365*8),50,cols(1,:),'filled'); text(pmlIII(doy+365*8)+0.5,zmlIII(doy+365*8),num2str(doy))
doy = 153;
scatter(pmlII(doy+365*8),zmlII(doy+365*8),50,cols(2,:),'filled'); text(pmlII(doy+365*8)+0.5,zmlII(doy+365*8),'153')
scatter(pmlIII(doy+365*8),zmlIII(doy+365*8),50,cols(1,:),'filled'); text(pmlIII(doy+365*8)+0.5,zmlIII(doy+365*8),'153')
doy = 300;
scatter(pmlII(doy+365*8),zmlII(doy+365*8),50,cols(2,:),'filled'); text(pmlII(doy+365*8)+0.5,zmlII(doy+365*8),'300')
scatter(pmlIII(doy+365*8),zmlIII(doy+365*8),50,cols(1,:),'filled'); text(pmlIII(doy+365*8)+0.5,zmlIII(doy+365*8),'300')
xlabel({'surface phytoplankton concentration (mg C m^{-3})'})
ylabel({'surface zooplankton concentration (mg C m^{-3})'})
set(gca,'fontsize',20)
hold off