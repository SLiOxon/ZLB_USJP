close all;
clear all;
clc;

%% US IRF (Figure 5)
Dinit = readmatrix("IRFlimits_US.csv"); %file must have no header and different series must be in columns
Dinit = Dinit(:, 2:end);

D(:,1:5) = Dinit(:,1:5);
D(:,6:9) = Dinit(:,14:17);
D(:,10:13) = Dinit(:,6:9);
D(:,14:17) = Dinit(:,18:21);
D(:,18:21) = Dinit(:,10:13);
D(:,22:25) = Dinit(:,22:25);

x = D(:,1)';

Dinit = readmatrix("IRFy2star_US.csv");
Dinit = Dinit(:, 2:end);
D(:,26:33) = Dinit(:,2:9);

% D_lim 8 by 2, first column min, second column max, first 4 rows for 1999,
% second 4 rows for 2009
D_lim = [min(D(:, [4 8 12 16 20 24 28 32])', [], 2), max(D(:, [5 9 13 17 21 25 29 33])', [], 2)];							% store limits for next plot
D_lim(:,1) = (D_lim(:,1) >= 0).*D_lim(:,1)/1.1 + (D_lim(:,1) < 0).*D_lim(:,1)*1.1;
D_lim(:,2) = (D_lim(:,2) >= 0).*D_lim(:,2)*1.1 + (D_lim(:,2) < 0).*D_lim(:,2)/1.1;
D_lim([5 6 7 8], 1) = -0.25;

areacol = [0 0 1]/1.4;
areacol_s = [1 1 1];
transp = 0.3;
linwid = 3;
axiswid = 1;

figure

for i = 1:8
subplot(4,2,i);
plot(x,0*x,'k','LineWidth',axiswid);
hold on
le = plot(x,D(:,4*i)',':r','LineWidth',axiswid,'DisplayName','67% error bands');
plot(x,D(:,4*i+1)',':r','LineWidth',axiswid);

xlim([min(x) max(x)]);
i_lim = -(floor(i/2)==ceil(i/2)) + i;
y_lo = min(D_lim(i_lim,1), D_lim(i_lim+1,1));
y_hi = max(D_lim(i_lim,2), D_lim(i_lim+1,2));
ylim([y_lo y_hi]);

p=patch([x fliplr(x)], [D(:,4*i-2)' fliplr(D(:,4*i-1)')], areacol);
p.FaceAlpha = transp;
p.EdgeAlpha = 0;

if i == 1 
    title('IRF to monetary policy shock in 1999Q1'); 
    ylabel('Inflation');
    legend(le);
end
if i == 2 
    title('IRF to monetary policy shock in 2009Q1'); 
    legend(le);
end
if i == 3 
    ylabel('Output Gap'); 
end
if i == 5 
    ylabel('Federal Funds Rate'); 
end
if i == 7 
    ylabel('Shadow Rate'); 
end

end

%% Japan IRF (Figure 6)

Dinit = readmatrix("IRFlimits_JP.csv"); %file must have no header and different series must be in columns
Dinit = Dinit(:, 2:end);

D(:,1:5) = Dinit(:,1:5);
D(:,6:9) = Dinit(:,14:17);
D(:,10:13) = Dinit(:,6:9);
D(:,14:17) = Dinit(:,18:21);
D(:,18:21) = Dinit(:,10:13);
D(:,22:25) = Dinit(:,22:25);

x = D(:,1)';

Dinit = readmatrix("IRFy2star_JP.csv");
Dinit = Dinit(:, 2:end);
D(:,26:33) = Dinit(:,2:9);

% D_lim 6 by 2, first column min, second column max, first 3 rows for 1999,
% second 3 rows for 2009
D_lim = [min(D(:, [4 8 12 16 20 24 28 32])', [], 2), max(D(:, [5 9 13 17 21 25 29 33])', [], 2)];							% store limits for next plot
D_lim(:,1) = (D_lim(:,1) >= 0).*D_lim(:,1)/1.1 + (D_lim(:,1) < 0).*D_lim(:,1)*1.1;
D_lim(:,2) = (D_lim(:,2) >= 0).*D_lim(:,2)*1.1 + (D_lim(:,2) < 0).*D_lim(:,2)/1.1;
D_lim([5 6 7 8], 1) = -0.25;

areacol = [0 0 1]/1.4;
areacol_s = [1 1 1];
transp = 0.3;
linwid = 3;
axiswid = 1;

figure

for i = 1:8
subplot(4,2,i);
plot(x,0*x,'k','LineWidth',axiswid);
hold on
le = plot(x,D(:,4*i)',':r','LineWidth',axiswid, 'DisplayName','67% error bands');
plot(x,D(:,4*i+1)',':r','LineWidth',axiswid);

xlim([min(x) max(x)]);
i_lim = -(floor(i/2)==ceil(i/2)) + i;
y_lo = min(D_lim(i_lim,1), D_lim(i_lim+1,1));
y_hi = max(D_lim(i_lim,2), D_lim(i_lim+1,2));
ylim([y_lo y_hi]);

p=patch([x fliplr(x)], [D(:,4*i-2)' fliplr(D(:,4*i-1)')], areacol);
p.FaceAlpha = transp;
p.EdgeAlpha = 0;

if i == 1 
    title('IRF to monetary policy shock in 1990Q1'); 
    ylabel('Inflation');
    legend(le);
end
if i == 2 
    title('IRF to monetary policy shock in 2010Q1'); 
    legend(le);
end
if i == 3 
    ylabel('Output Gap'); 
end
if i == 5 
    ylabel('Call Rate'); 
end
if i == 7 
    ylabel('Shadow Rate'); 
end

end

%% US Cumulative (Figure 7)

D = readmatrix("IRFlimits_horz_US.csv"); %file must have no header and different series must be in columns
D = D(:, 2:end);
D(131:135, [3:4,9:10,15:16,21:22]) = NaN;
D(127:135, [5:6,11:12,17:18,23:24]) = NaN;

% D_lim 9 by 2, first column min, second column max, first 3 rows for 1999,
% second 3 rows for 2009
D_lim = [min(D(:, [1 3 5 7 9 11 13 15 17 19 21 23])', [], 2), max(D(:, [2 4 6 8 10 12 14 16 18 20 22 24])', [], 2)];							% store limits for next plot
D_lim(:,1) = (D_lim(:,1) >= 0).*D_lim(:,1)/1.1 + (D_lim(:,1) < 0).*D_lim(:,1)*1.1;
D_lim(:,2) = (D_lim(:,2) >= 0).*D_lim(:,2)*1.1 + (D_lim(:,2) < 0).*D_lim(:,2)/1.1;

areacol = [0 0 1]/1.4;
areacol_s = [1 1 1];
transp = 0.3;
linwid = 3;
axiswid = 1;

figure
starttime = datenum('01-jul-1985','dd-mmm-yyyy');
endtime = datenum('01-jan-2019','dd-mmm-yyyy');
xdata = linspace(starttime, endtime, 135);
endtime2 = datenum('01-oct-2017','dd-mmm-yyyy');
xdata2 = linspace(starttime, endtime, 130);
endtime3 = datenum('01-oct-2016','dd-mmm-yyyy');
xdata3 = linspace(starttime, endtime, 126);

for i = 1:12
subplot(4,3,i);

if (i==1)|(i==4)|(i==7)|(i==10)
    plot(xdata,0*xdata,'k','LineWidth',axiswid);
    hold on

    xlim([min(xdata) max(xdata)]);
    ylim([D_lim(i,1) D_lim(i,2)]);
    
    p=patch([xdata fliplr(xdata)], [D(:,2*i-1)' fliplr(D(:,2*i)')], areacol);
    p.FaceAlpha = transp;
    p.EdgeAlpha = 0;
    
    xticks([datenum('01-jan-1990','dd-mmm-yyyy') datenum('01-jan-2000','dd-mmm-yyyy') datenum('01-jan-2010','dd-mmm-yyyy')]);
    datetick('x','yyyy','keeplimits','keepticks');
end
if (i==2)|(i==5)|(i==8)|(i==11)
    plot(xdata2,0*xdata2,'k','LineWidth',axiswid);
    hold on

    xlim([min(xdata2) max(xdata2)]);
    ylim([D_lim(i,1) D_lim(i,2)]);
    
    p=patch([xdata2 fliplr(xdata2)], [D(1:130,2*i-1)' fliplr(D(1:130,2*i)')], areacol);
    p.FaceAlpha = transp;
    p.EdgeAlpha = 0;
    
    xticks([datenum('01-jan-1990','dd-mmm-yyyy') datenum('01-jan-2000','dd-mmm-yyyy') datenum('01-jan-2010','dd-mmm-yyyy')]);
    datetick('x','yyyy','keeplimits','keepticks');
end
if (i==3)|(i==6)|(i==9)|(i==12)
    plot(xdata3,0*xdata3,'k','LineWidth',axiswid);
    hold on

    xlim([min(xdata3) max(xdata3)]);
    ylim([D_lim(i,1) D_lim(i,2)]);
    
    p=patch([xdata3 fliplr(xdata3)], [D(1:126,2*i-1)' fliplr(D(1:126,2*i)')], areacol);
    p.FaceAlpha = transp;
    p.EdgeAlpha = 0;
    
    xticks([datenum('01-jan-1990','dd-mmm-yyyy') datenum('01-jan-2000','dd-mmm-yyyy') datenum('01-jan-2010','dd-mmm-yyyy')]);
    datetick('x','yyyy','keeplimits','keepticks');
end
% p.FaceAlpha = transp;
% p.EdgeAlpha = 0;
% candystripe(p, 'Color',areacol_s,'Width',1);

if i == 1 
    title('On impact'); 
    ylabel('Inflation');
end
if i == 2 
    title('4 quarters ahead'); 
end
if i == 3 
    title('8 quarters ahead'); 
end
if i == 4
    ylabel('Output Gap');
end
if i == 7
    ylabel('Federal Funds Rate'); 
end
if i == 10
    ylabel('Shadow Rate'); 
end

end

%% Japan Cumulative (Figure 8)

D = readmatrix("IRFlimits_horz_JP.csv"); %file must have no header and different series must be in columns
D = D(:, 2:end);
D(131:135, [3:4,9:10,15:16,21:22]) = NaN;
D(127:135, [5:6,11:12,17:18,23:24]) = NaN;

% D_lim 9 by 2, first column min, second column max, first 3 rows for 1999,
% second 3 rows for 2009
D_lim = [min(D(:, [1 3 5 7 9 11 13 15 17 19 21 23])', [], 2), max(D(:, [2 4 6 8 10 12 14 16 18 20 22 24])', [], 2)];					% store limits for next plot
D_lim(:,1) = (D_lim(:,1) >= 0).*D_lim(:,1)/1.1 + (D_lim(:,1) < 0).*D_lim(:,1)*1.1;
D_lim(:,2) = (D_lim(:,2) >= 0).*D_lim(:,2)*1.1 + (D_lim(:,2) < 0).*D_lim(:,2)/1.1;

areacol = [0 0 1]/1.4;
areacol_s = [1 1 1];
transp = 0.3;
linwid = 3;
axiswid = 1;

figure
starttime = datenum('01-jul-1985','dd-mmm-yyyy');
endtime = datenum('01-jan-2019','dd-mmm-yyyy');
xdata = linspace(starttime, endtime, 135);
endtime2 = datenum('01-oct-2017','dd-mmm-yyyy');
xdata2 = linspace(starttime, endtime, 130);
endtime3 = datenum('01-oct-2016','dd-mmm-yyyy');
xdata3 = linspace(starttime, endtime, 126);

for i = 1:12
subplot(4,3,i);

if (i==1)|(i==4)|(i==7)|(i==10)
    plot(xdata,0*xdata,'k','LineWidth',axiswid);
    hold on

    xlim([min(xdata) max(xdata)]);
    ylim([D_lim(i,1) D_lim(i,2)]);
    
    p=patch([xdata fliplr(xdata)], [D(:,2*i-1)' fliplr(D(:,2*i)')], areacol);
    p.FaceAlpha = transp;
    p.EdgeAlpha = 0;
    
    xticks([datenum('01-jan-1990','dd-mmm-yyyy') datenum('01-jan-2000','dd-mmm-yyyy') datenum('01-jan-2010','dd-mmm-yyyy')]);
    datetick('x','yyyy','keeplimits','keepticks');
end
if (i==2)|(i==5)|(i==8)|(i==11)
    plot(xdata2,0*xdata2,'k','LineWidth',axiswid);
    hold on

    xlim([min(xdata2) max(xdata2)]);
    ylim([D_lim(i,1) D_lim(i,2)]);
    
    p=patch([xdata2 fliplr(xdata2)], [D(1:130,2*i-1)' fliplr(D(1:130,2*i)')], areacol);
    p.FaceAlpha = transp;
    p.EdgeAlpha = 0;
    
    xticks([datenum('01-jan-1990','dd-mmm-yyyy') datenum('01-jan-2000','dd-mmm-yyyy') datenum('01-jan-2010','dd-mmm-yyyy')]);
    datetick('x','yyyy','keeplimits','keepticks');
end
if (i==3)|(i==6)|(i==9)|(i==12)
    plot(xdata3,0*xdata3,'k','LineWidth',axiswid);
    hold on

    xlim([min(xdata3) max(xdata3)]);
    ylim([D_lim(i,1) D_lim(i,2)]);
    
    p=patch([xdata3 fliplr(xdata3)], [D(1:126,2*i-1)' fliplr(D(1:126,2*i)')], areacol);
    p.FaceAlpha = transp;
    p.EdgeAlpha = 0;
    
    xticks([datenum('01-jan-1990','dd-mmm-yyyy') datenum('01-jan-2000','dd-mmm-yyyy') datenum('01-jan-2010','dd-mmm-yyyy')]);
    datetick('x','yyyy','keeplimits','keepticks');
end
% p.FaceAlpha = transp;
% p.EdgeAlpha = 0;
% candystripe(p, 'Color',areacol_s,'Width',1);

if i == 1 
    title('On impact'); 
    ylabel('Inflation');
end
if i == 2 
    title('4 quarters ahead'); 
end
if i == 3 
    title('8 quarters ahead'); 
end
if i == 4
    ylabel('Output Gap');
end
if i == 7
    ylabel('Call Rate'); 
end
if i == 10
    ylabel('Shadow Rate'); 
end

end

%% US Shadow Rate (Figure 9)

D = readmatrix("shadowlimits_US.csv"); %file must have no header and different series must be in columns
D = D(:, 2:end);

areacol = [0 0 1]/1.4;
areacol_s = [1 1 1];
transp = 0.3;
linwid = 3;
axiswid = 1;

figure
starttime = datenum('01-jul-1985','dd-mmm-yyyy');
endtime = datenum('01-jan-2019','dd-mmm-yyyy');
xdata = linspace(starttime, endtime, 135);

plot(xdata,0*xdata,'k','LineWidth',axiswid);
hold on

xlim([min(xdata) max(xdata)]);
    
p=patch([xdata fliplr(xdata)], [D(:,1)' fliplr(D(:,2)')], areacol);
p.FaceAlpha = transp;
p.EdgeAlpha = 1;
    
datetick('x','yyyy','keeplimits');

%% Japan Shadow Rate (Figure 10)

D = readmatrix("shadowlimits_JP.csv"); %file must have no header and different series must be in columns
D = D(:, 2:end);

areacol = [0 0 1]/1.4;
areacol_s = [1 1 1];
transp = 0.3;
linwid = 3;
axiswid = 1;

figure
starttime = datenum('01-jul-1985','dd-mmm-yyyy');
endtime = datenum('01-jan-2019','dd-mmm-yyyy');
xdata = linspace(starttime, endtime, 135);

plot(xdata,0*xdata,'k','LineWidth',axiswid);
hold on

xlim([min(xdata) max(xdata)]);
    
p=patch([xdata fliplr(xdata)], [D(:,1)' fliplr(D(:,2)')], areacol);
p.FaceAlpha = transp;
p.EdgeAlpha = 1;
    
datetick('x','yyyy','keeplimits');