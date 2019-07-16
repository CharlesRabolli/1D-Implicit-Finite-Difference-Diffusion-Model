%% Bolus Injection vs Construct

clear all; %Deletes anything previously saved
clc; %clears the command window
InitialDrugConcentration = 1; %mM, initial concentration that we are using
InitialSA = 0.332; % cm^2, surcface area of our plate
InitialHeight = .13; %cm, height of our model
InitialVolume = InitialSA*InitialHeight; %cm^3, volume of model
InitialMolesDrug = 1*10^12; %femtomoles
D_Media = 0.024156; %cm^2/hr, found online
L1 = 0.04; %cm - distance between transwell to cells
T1 = 96; %hour, arbitrary time 
dx1=(L1)/10; % Step size in the x direction
tstepsize = 10000; % number of points
dt1=1/tstepsize; % Step size in the t direction
xmesh1=0:dx1:L1; % Creates an array of values in the x dir.
tmesh1=0:dt1:T1; % Creates an array of values in the t dir.
nx1 = length(xmesh1); % number of points in x dimension
nt1 = length(tmesh1); % number of points in t dimension
detfunc = (nt1-1)/(tstepsize*8); %We used 8 hour time points;
lambda1 = D_Media * dt1 / (dx1^2); % stepsize for numerical integration
totalfunc = zeros([1,length(tmesh1)]); % create matrix of zeroes
Interface = zeros([1,length(tmesh1)]); % create matrix of zeroes
Intsum = zeros([1,length(tmesh1)]); % create matrix of zeroes 
solution1 = zeros(nt1, nx1); % Creates a matrix size x by t to find the diffusion before beads
AmountReleased = [0, 0.1, 0.072, 0.0752, 0.05, 0.062, 0.0567, 0.0743, 0.0757, 0.0711, 0.06055 0.05, 0.01].*InitialMolesDrug; %these numbers were determined experimentally
TimelyRelease = AmountReleased/((nt1-1)/(length(AmountReleased)-1)); % creates a vector to distribute the amount over time
for n =1:T1*tstepsize %Moves forward in the t direction until the last time point
StepTimelyRelease(n) = TimelyRelease(ceil(n/80000)+1);
end
for t1 = 2:nt1
    solution1(t1,1)=solution1(t1,1)+ StepTimelyRelease(t1-1);
    solution1(t1,nx1)=solution1(t1-1,nx1-1);
 for x1 = 2:nx1-1
  solution1(t1,x1) = solution1(t1-1,x1) + lambda1 * ...
 (solution1(t1-1,x1+1) - 2 * solution1(t1-1, x1) + solution1(t1-1,x1-1));
     Interface(t1) = solution1(t1,1);
     Intsum(t1) = Interface(t1) + Intsum(t1-1);
 end
solution1(t1,2) = solution1(t1, 2)*exp(-0.256721*(rem(t1-1,tstepsize)/tstepsize));
     totalfunc(t1) = solution1(t1,nx1-1)+totalfunc(t1-1);
end
v = size(totalfunc);% 1 row x ? columns, but we only care about columns so...
qwe = v(2); % Now we know how many columns there are
qwestep = (qwe-1)/tstepsize; %Matlab indexes at 0, so subtract one
qweint = 1:1/tstepsize:qwestep+1; % creates a vector until the end, by our interval, within our time range.


lambda2 = D_Media * dt1 / (dx1^2); % stepsize for numerical integration
totalfunc2 = zeros([1,length(tmesh1)]);
Interface2 = zeros([1,length(tmesh1)]);
Intsum2 = zeros([1,length(tmesh1)]);
solution2 = zeros(nt1, nx1); % Creates a matrix size x by t to find the diffusion before beads
for t1 = 2:nt1
    solution2(2,1)=InitialMolesDrug;
    solution2(t1,nx1)=solution2(t1-1,nx1-1);
 for x1 = 2:nx1-1
  solution2(t1,x1) = solution2(t1-1,x1) + lambda2 * ...
 (solution2(t1-1,x1+1) - 2 * solution2(t1-1, x1) + solution2(t1-1,x1-1));
     Interface2(t1) = solution2(t1,1);
     Intsum2(t1) = Interface2(t1) + Intsum2(t1-1);
 end
solution2(t1,2) = solution2(t1, 2)*exp(-0.256721*(rem(t1-1,tstepsize)/tstepsize));
     totalfunc2(t1) = solution2(t1,nx1-1)+totalfunc2(t1-1);
end
figure(1)
v2 = size(totalfunc2);
qwe2 = v(2);
qwestep2 = (qwe-1)/tstepsize;
qweint2 = 1:1/tstepsize:qwestep+1;
plot(qweint,totalfunc/InitialMolesDrug*100)
hold on 
plot(qweint2,totalfunc2/InitialMolesDrug*100)
ylabel('Percent')
xlabel('Hours')
title('Percent of Initial 1 millimole that Reaches Beads')
axis([0 100 0 100])

HourSum1 = zeros([1,96]);
HourSum2 = zeros([1,96]);
k = 1;
span = 1:96;
hourly = span * 10000;
dist = 1:960000;
for j = 1:nt1-1
    if j < hourly(ceil(j/10000))
        HourSum1(k) = solution1(j,nx1-1) + HourSum1(k);
        HourSum2(k) = solution2(j,nx1-1) + HourSum2(k);
    elseif j == hourly(ceil(j/10000))
        k = k+1;
    end
end

span2 = 8:8:96;
THourSum1(1) = sum(HourSum1(1:8)); % Create a new number, that is the sum of everything in that time range.
THourSum1(2) = sum(HourSum1(9:16));
THourSum1(3) = sum(HourSum1(17:24));
THourSum1(4) = sum(HourSum1(25:32));
THourSum1(5) = sum(HourSum1(33:40));
THourSum1(6) = sum(HourSum1(41:48));
THourSum1(7) = sum(HourSum1(49:56));
THourSum1(8) = sum(HourSum1(57:64));
THourSum1(9) = sum(HourSum1(65:72));
THourSum1(10) = sum(HourSum1(73:80));
THourSum1(11) = sum(HourSum1(81:88));
THourSum1(12) = sum(HourSum1(89:96));
THourSum2(1) = sum(HourSum2(1:8));
THourSum2(2) = sum(HourSum2(9:16));
THourSum2(3) = sum(HourSum2(17:24));
THourSum2(4) = sum(HourSum2(25:32));
THourSum2(5) = sum(HourSum2(33:40));
THourSum2(6) = sum(HourSum2(41:48));
THourSum2(7) = sum(HourSum2(49:56));
THourSum2(8) = sum(HourSum2(57:64));
THourSum2(9) = sum(HourSum2(65:72));
THourSum2(10) = sum(HourSum2(73:80));
THourSum2(11) = sum(HourSum2(81:88));
THourSum2(12) = sum(HourSum2(89:96));
figure(2)
hold on % plot figures on top of one another
plot(span2,THourSum1/InitialMolesDrug*100) % create a percent
plot(span2,THourSum2/InitialMolesDrug*100)
set(gca,'XTick',(0:8:96)) % Change axis interval
title('Amount of Drug Reaching Beads Every 8 Hours')
xlabel('Hours')
ylabel('Percent of Initial Dose')

InitialMolesDrugCons = 1*10^(12) * 0.430012; %femtomoles
D_MediaAlg = 9.972E-4; %cm^2/hr
L1Cell = 0.0225; %cm - distance between transwell to cells
T1 = 96; %hour
dx1=(L1Cell)/10; % Step size in the x direction
tstepsize = 10000; %previously 10000
dt1=1/tstepsize; % Step size in the t direction
xmesh1=0:dx1:L1Cell; % Creates an array of values in the x dir.
tmesh1=0:dt1:T1; % Creates an array of values in the t dir.
nx1 = length(xmesh1); % number of points in x dimension
nt1 = length(tmesh1); % number of points in t dimension
lambda1Cell = D_MediaAlg * (dt1 / ((dx1)^2)); % stepsize for numerical integration
totalfunc1Cell = zeros([1,length(tmesh1)]);
Interface1Cell = zeros([1,length(tmesh1)]);
Intsum1Cell = zeros([1,length(tmesh1)]);
solution1Cell = zeros(nt1, nx1); % Creates a matrix size x by t to find the diffusion before beads
AmountReleased1Cell = [0,   5.6628    4.0906    4.2681    2.8415    3.5176    3.2191    4.2150    4.2967    4.0364    3.4383    2.8395    0.5730 ].*10^(10);
TimelyRelease1Cell = AmountReleased1Cell/((nt1-1)/(length(AmountReleased1Cell)-1));
for n =1:T1*tstepsize
StepTimelyRelease1Cell(n) = TimelyRelease1Cell(ceil(n/80000)+1);%previously 80000
end
TotalStepRelease1Cell = zeros([1,length(tmesh1)]);
for t1 = 2:nt1
    solution1Cell(t1,1)=StepTimelyRelease1Cell(t1-1);
    solution1Cell(t1,nx1)=solution1Cell(t1-1,nx1-1);
 for x1 = 2:nx1-1
  solution1Cell(t1,x1) = solution1Cell(t1-1,x1) + lambda1Cell * ...
 (solution1Cell(t1-1,x1+1) - 2 * solution1Cell(t1-1, x1) + solution1Cell(t1-1,x1-1));
     Interface1Cell(t1) = solution1Cell(t1,1);
     Intsum1Cell(t1) = Interface1Cell(t1) + Intsum1Cell(t1-1);
 end
solution1Cell(t1,2) = solution1Cell(t1, 2)*exp(-0.256721*(rem(t1-1,tstepsize)/tstepsize));
   totalfunc1Cell(t1) = solution1Cell(t1,nx1-1)+totalfunc1Cell(t1-1);
end
figure
v1Cell = size(totalfunc1Cell);
qwe1Cell = v1Cell(2);
qwestep1Cell = (qwe1Cell-1)/tstepsize;
qweint1Cell = 1:1/tstepsize:qwestep1Cell+1;
plot(qweint1Cell,totalfunc1Cell/(10^12)*100)
ylabel('Percent')
xlabel('Hours')
title('Percent of Initial 1 millimole that Reaches Beads')
axis([0 100 0 10])
ReachCells1Cell = totalfunc1Cell(nt1);
ReachCells1Cell
PercentToCells1Cell = ReachCells1Cell/(10^12) * 100
TotalBoundary11Cell = sum(Interface1Cell)
Boundary1Cell = Interface1Cell/InitialMolesDrugCons * 100;
figure
plot(qweint1Cell,totalfunc1Cell/InitialMolesDrugCons*100)
ylabel('Percent')
xlabel('Hours')
title('Percent of Drug Entering Beads that reaches cells')
axis([0 100 0 25])

InitialMolesDrugBol = 1*10^12 * 0.990238; %femtomoles
D_Alginate = 9.972E-4; %cm^2/hr
L2Cell = 0.0225; %cm - distance between transwell to cells
T1 = 96; %hour
dx1=(L2Cell)/10; % Step size in the x direction
tstepsize = 10000;
dt1=1/tstepsize; % Step size in the t direction
xmesh1=0:dx1:L2Cell; % Creates an array of values in the x dir.
tmesh1=0:dt1:T1; % Creates an array of values in the t dir.
nx1 = length(xmesh1); % number of points in x dimension
nt1 = length(tmesh1); % number of points in t dimension
lambda2Cell = D_Alginate * dt1 / (dx1^2); % stepsize for numerical integration
totalfunc2Cell = zeros([1,length(tmesh1)]);
Interface2Cell = zeros([1,length(tmesh1)]);
Intsum2Cell = zeros([1,length(tmesh1)]);
solution2Cell = zeros(nt1, nx1); % Creates a matrix size x by t to find the diffusion before beads
AmountReleased2Cell = [0, 9.90238, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0].*(1*10^11);
TimelyRelease2Cell = AmountReleased2Cell/((nt1-1)/(length(AmountReleased2Cell)-1));
for n =1:T1*tstepsize
StepTimelyRelease2Cell(n) = TimelyRelease2Cell(ceil(n/80000)+1);
end
for t1 = 2:nt1
       solution2Cell(t1,1)=StepTimelyRelease2Cell(t1-1);
    solution2Cell(t1,nx1)=solution2Cell(t1-1,nx1-1);
 for x1 = 2:nx1-1
  solution2Cell(t1,x1) = solution2Cell(t1-1,x1) + lambda2Cell * ...
 (solution2Cell(t1-1,x1+1) - 2 * solution2Cell(t1-1, x1) + solution2Cell(t1-1,x1-1));
     Interface2Cell(t1) = solution2Cell(t1,1);
     Intsum2Cell(t1) = Interface2Cell(t1) + Intsum2Cell(t1-1);
 end
solution2Cell(t1,2) = solution2Cell(t1, 2)*exp(-0.256721*(rem(t1-1,tstepsize)/tstepsize));
     totalfunc2Cell(t1) = solution2Cell(t1,nx1-1)+totalfunc2Cell(t1-1);
end
figure
v2Cell = size(totalfunc2Cell);
qwe2Cell = v2Cell(2);
qwestep2Cell = (qwe2Cell-1)/tstepsize;
qweint2Cell = 1:1/tstepsize:qwestep2Cell+1;
plot(qweint2Cell,totalfunc2Cell/InitialMolesDrugBol*100)
ylabel('Percent')
xlabel('Hours')
title('Percent of Drug Entering Beads that reaches cells')
axis([0 100 0 100])
ReachCells2Cell = totalfunc2Cell(nt1);
ReachCells2Cell
PercentToCells2Cell = ReachCells2Cell/(1*10^12) * 100
TotalBoundary2Cell = sum(Interface2Cell)
Boundary2Cell = Interface2Cell/(1*10^12) * 100;
figure
plot(qweint2Cell,totalfunc2Cell/(1*10^12)*100, 'r')
ylabel('Percent')
xlabel('Hours')
title('Percent of Initial 1mmole that reaches cells')
axis([0 100 0 20])

HourSum1Cell = zeros([1,96]);
HourSum2Cell = zeros([1,96]);
k = 1;
span = 1:96;
hourly = span * 10000;
dist = 1:960000;
for j = 1:nt1-1
    if j < hourly(ceil(j/10000))
        HourSum1Cell(k) = solution1Cell(j,nx1-1) + HourSum1Cell(k);
        HourSum2Cell(k) = solution2Cell(j,nx1-1) + HourSum2Cell(k);
    elseif j == hourly(ceil(j/10000))
        k = k+1;
    end
end
span2 = 8:8:96;
THourSum1Cell(1) = sum(HourSum1Cell(1:8));
THourSum1Cell(2) = sum(HourSum1Cell(9:16));
THourSum1Cell(3) = sum(HourSum1Cell(17:24));
THourSum1Cell(4) = sum(HourSum1Cell(25:32));
THourSum1Cell(5) = sum(HourSum1Cell(33:40));
THourSum1Cell(6) = sum(HourSum1Cell(41:48));
THourSum1Cell(7) = sum(HourSum1Cell(49:56));
THourSum1Cell(8) = sum(HourSum1Cell(57:64));
THourSum1Cell(9) = sum(HourSum1Cell(65:72));
THourSum1Cell(10) = sum(HourSum1Cell(73:80));
THourSum1Cell(11) = sum(HourSum1Cell(81:88));
THourSum1Cell(12) = sum(HourSum1Cell(89:96));
THourSum2Cell(1) = sum(HourSum2Cell(1:8));
THourSum2Cell(2) = sum(HourSum2Cell(9:16));
THourSum2Cell(3) = sum(HourSum2Cell(17:24));
THourSum2Cell(4) = sum(HourSum2Cell(25:32));
THourSum2Cell(5) = sum(HourSum2Cell(33:40));
THourSum2Cell(6) = sum(HourSum2Cell(41:48));
THourSum2Cell(7) = sum(HourSum2Cell(49:56));
THourSum2Cell(8) = sum(HourSum2Cell(57:64));
THourSum2Cell(9) = sum(HourSum2Cell(65:72));
THourSum2Cell(10) = sum(HourSum2Cell(73:80));
THourSum2Cell(11) = sum(HourSum2Cell(81:88));
THourSum2Cell(12) = sum(HourSum2Cell(89:96));
figure
hold on
plot(span2,THourSum1Cell/InitialMolesDrug*100)
plot(span2,THourSum2Cell/InitialMolesDrug*100)
set(gca,'XTick',(0:8:96))
title('Amount of Drug Reaching Cells Every 8 Hours')
xlabel('Hours')
ylabel('Percent of Initial Dose')
axis([0 100 0 20])

%% Using BreakPlot - Be careful! This alters the data point! Need to run original code immediately before running this
figure
hPlotData = breakplot(span2,THourSum2Cell/InitialMolesDrug*100,2,16,'',2);
delete(hPlotData);
THourSum2Cell(1) = 0.4460* 10^11;
hold on
plot(span2,THourSum1Cell/InitialMolesDrug*100,'color','k','linewidth',1.5)
plot(span2,THourSum2Cell/InitialMolesDrug*100,'color','k', 'linestyle','-.','linewidth',1)
plot(span2,THourSum1Cell/InitialMolesDrug*100,'^','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5)
plot(span2,THourSum2Cell/InitialMolesDrug*100,'d','MarkerFaceColor','k','MarkerEdgeColor','k','MarkerSize',5)
set(gca,'XTick',(0:8:96))
set(gca, 'fontsize',11)
set(gca, 'Fontname','times')
set(gca, 'box', 'off')
xlabel('Time (Hr)')
ylabel('Percent of Initial Dose')
%axis([0 100 0 50])
legend('Construct','Bolus','location','NW')
legend boxoff