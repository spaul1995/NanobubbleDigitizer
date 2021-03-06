% NanobubbleDigitizer data analysis package licence details:
% 
% Copyright 2009-2014 Soumyadeep Paul (soumyadeep.paul@thml.t.u-tokyo.ac.jp)
% 
% This file is part of NanobubbleDigitizer data analysis package.
% 
% NanobubbleDigitizer is a MATLAB package designed to detect discrete 
% bubble events from transient current traces. The bubble features
% extracted from the current data using this code can be used to analyze
% nanoscopic bubble dynamics confined within solid-state nanopores. This
% package is released under the GNU GPL.
% 
% The copyright of NanobubbleDigitizer is held by Soumyadeep Paul.
% 
% NanobubbleDigitizer is released under the GNU GPL.  NanobubbleDigitizer 
% is free software: you can redistribute it and/or modify it under the 
% terms of the GNU General Public License (version 3) as published 
% by the Free Software Foundation.You should have received a copy of the 
% GNU General Public Licence along with arb (see file licence.txt after 
% unpacking).  If not, see <http://www.gnu.org/licences/>.
% 
% NanobubbleDigitizer is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General 
% Public Licence for more details.
% 
% For full details of NanobubbleDigitizer's licence see the licence file in the main 
% directory.

clear all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section: 0   Control parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counterpulse = 1;
filename = sprintf('M%i.csv',fix(counterpulse));
start = 0;
%ending = 137600;
ending = 246250;
%ending = 621150;
X1 = csvread(filename,start,0,[start,0,ending,3]);

startpt = 200*(counterpulse-1);
endptt = 100+200*(counterpulse-1);
T1 = 0.0004;
RawSignal = X1(:,2);
LowFilteredSignal = X1(:,3);
t = X1(:,1);

base1 = 150;
base2 = 350;

base12 = 200;
base15 = 250;
base18 = 300;

K = 10e5;
s1 = 950;
s2 = s1-450;
cutCurrentd = @(t) s1-base1-s2*(1-exp(-K*1e-6*(t-startpt)));
cutCurrent12 = @(t) s1-base12-s2*(1-exp(-K*1e-6*(t-startpt)));
cutCurrent15 = @(t) s1-base15-s2*(1-exp(-K*1e-6*(t-startpt)));
cutCurrent18 = @(t) s1-base18-s2*(1-exp(-K*1e-6*(t-startpt)));
cutCurrentr = @(t) s1-base2-s2*(1-exp(-K*1e-6*(t-startpt)));


sigma = 5;
recede_fac_ini = 0.01;
order = 3;
averagingfac=20/1250;
nucpoint_averaging_width = 28;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Section: 1   Start of Cut point estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cutpt1 = zeros(1000,1);
Cutpt12 = zeros(1000,1);
Cutpt15 = zeros(1000,1);
Cutpt18 = zeros(1000,1);
Cutpt2dash = zeros(1000,1);

counter1 = 1;
counter12 = 1;
counter15 = 1;
counter18 = 1;
counter2 = 1;
i = 1;
spacer = 10;
switch1 = 0;
while i < length(RawSignal)
   if (RawSignal(i)<cutCurrentd(t(i))+sigma/2) && (RawSignal(i)>cutCurrentd(t(i))-sigma/2) && sign(RawSignal(i+200)-RawSignal(i)) < 0 && sign(RawSignal(i+50)-RawSignal(i)) < 0 && switch1 == 0
       Cutpt1(counter1) = i;
       counter1 = counter1+1;
       i = i+spacer;
       switch2 = 0;
       switch1 = 1;
       switch12 = 0;
       switch15 = 0;
       switch18 = 0;
       
   end
   
   
   if (RawSignal(i)<cutCurrent12(t(i))+sigma/2) && (RawSignal(i)>cutCurrent12(t(i))-sigma/2) && sign(RawSignal(i+spacer)-RawSignal(i)) > 0 && switch12 == 0
       Cutpt12(counter1) = i;
       counter12 = counter12+1;
       i = i+spacer;
       switch2 = 0;
       switch1 = 0;
       switch12 = 1;
       switch15 = 0;
       switch18 = 0;
   end
   if (RawSignal(i)<cutCurrent15(t(i))+sigma/2) && (RawSignal(i)>cutCurrent15(t(i))-sigma/2) && sign(RawSignal(i+spacer)-RawSignal(i)) > 0 && switch15 == 0
       Cutpt15(counter1) = i;
       counter15 = counter15+1;
       i = i+spacer;
       switch2 = 0;
       switch1 = 0;
       switch12 = 0;
       switch15 = 1;
       switch18 = 0;
   end
   if (RawSignal(i)<cutCurrent18(t(i))+sigma/2) && (RawSignal(i)>cutCurrent18(t(i))-sigma/2) && sign(RawSignal(i+spacer)-RawSignal(i)) > 0 && switch18 == 0
       Cutpt18(counter1) = i;
       counter18 = counter18+1;
       i = i+spacer;
       switch2 = 0;
       switch1 = 0;
       switch12 = 0;
       switch15 = 0;
       switch18 = 1;
   end
   
   
  if (RawSignal(i)<cutCurrentr(t(i))+sigma/2) && (RawSignal(i)>cutCurrentr(t(i))-sigma/2) && sign(RawSignal(i+spacer)-RawSignal(i)) > 0 && switch2 == 0
       Cutpt2dash(counter2) = i;
       counter2 = counter2+1;
       i = i+spacer;
       switch2 = 1;
       switch1 = 0;
       switch12 = 0;
       switch15 = 0;
       switch18 = 0;
  end
   i = i+1;
end

Cutpt2dash(Cutpt2dash==0)=[];
Cutpt1(Cutpt1==0)=[];

Cutpt12(Cutpt12==0)=[];
Cutpt15(Cutpt15==0)=[];
Cutpt18(Cutpt18==0)=[];

% len_1 = sprintf('lenght of Cutpt1 = %i',fix(length(Cutpt1)));
% len_12 = sprintf('lenght of Cutpt12 = %i',fix(length(Cutpt12)));
% len_15 = sprintf('lenght of Cutpt15 = %i',fix(length(Cutpt15)));
% len_18 = sprintf('lenght of Cutpt18 = %i',fix(length(Cutpt18)));
% len_2 = sprintf('lenght of Cutpt2dash = %i',fix(length(Cutpt2dash)));
% display(len_1)
% display(len_12)
% display(len_15)
% display(len_18)
% display(len_2)


if length(Cutpt2dash) > length(Cutpt1)
    display('error')
end

if length(Cutpt2dash) ~= length(Cutpt1)
figure
hold on
plot(t,RawSignal)
hold on
fplot(cutCurrentd,'-r')
hold on
fplot(cutCurrentr,'-g')
hold on
fplot(cutCurrent12,'-c')
hold on
fplot(cutCurrent15,'-k')
hold on
fplot(cutCurrent18,'-y')
hold on
plot(t(Cutpt1),RawSignal(Cutpt1),'*r')
hold on
plot(t(Cutpt2dash),RawSignal(Cutpt2dash),'*g')
hold on
plot(t(Cutpt12),RawSignal(Cutpt12),'*c')
hold on
plot(t(Cutpt15),RawSignal(Cutpt15),'*k')
hold on
plot(t(Cutpt18),RawSignal(Cutpt18),'*y')

hold off
xlim([startpt endptt])
ylim([0 1000])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of Cut point estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Section: 2   Choosing the effective 2nd Cut point from Cutpt12, Cutpt15, Cutpt18 and
%Cuttpt2dash
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lengthfit = zeros(length(Cutpt1),1)+70;
Cutpt2dashdash = zeros(length(Cutpt1),1);

for i=1:length(Cutpt1)
    if (Cutpt2dash(i) == min([Cutpt12(i),Cutpt15(i),Cutpt18(i), Cutpt2dash(i)])) && (RawSignal(Cutpt2dash(i)) > min([RawSignal(Cutpt1(i):Cutpt2dash(i))])+50)
        Cutpt2dashdash(i) = Cutpt2dash(i);
    elseif Cutpt18(i) == min([Cutpt12(i),Cutpt15(i),Cutpt18(i)]) && RawSignal(Cutpt18(i)) > min([RawSignal(Cutpt1(i):Cutpt18(i))])+50
        Cutpt2dashdash(i) = Cutpt18(i);
        
        if (Cutpt2dash(i) ~= min([Cutpt12(i),Cutpt15(i),Cutpt18(i), Cutpt2dash(i)]))
            b = [Cutpt2dash(1:i-1);0;Cutpt2dash(i:end)];
            Cutpt2dash = b;
            display('cutpt2dash miss')
            display(i)
        end
        
    elseif Cutpt15(i) == min(Cutpt12(i),Cutpt15(i)) && RawSignal(Cutpt15(i)) > min([RawSignal(Cutpt1(i):Cutpt15(i))])+50
        Cutpt2dashdash(i) = Cutpt15(i);
        
        if (Cutpt2dash(i) ~= min([Cutpt12(i),Cutpt15(i),Cutpt18(i), Cutpt2dash(i)]))
        b = [Cutpt2dash(1:i-1);0;Cutpt2dash(i:end)];
        Cutpt2dash = b;
        display('cutpt2dash miss')
        display(i)
        end
        
        if Cutpt18(i) ~= min([Cutpt12(i),Cutpt15(i),Cutpt18(i)])
        b = [Cutpt18(1:i-1);0;Cutpt18(i:end)];
        Cutpt18 = b;
        display('cutpt18 miss')
        display(i)
        end

    else
        Cutpt2dashdash(i) = Cutpt12(i);
        [minC,minindex] = min([RawSignal(Cutpt1(i):Cutpt12(i))]);
        lengthfit(i) = ceil(0.5*(length([RawSignal(Cutpt1(i):Cutpt12(i))])-minindex));
        
        if (Cutpt2dash(i) ~= min([Cutpt12(i),Cutpt15(i),Cutpt18(i), Cutpt2dash(i)]))
        b = [Cutpt2dash(1:i-1);0;Cutpt2dash(i:end)];
        Cutpt2dash = b;
        display('cutpt2dash miss')
        display(i)
        end
        
        if Cutpt18(i) ~= min([Cutpt12(i),Cutpt15(i),Cutpt18(i)])
        b = [Cutpt18(1:i-1);0;Cutpt18(i:end)];
        Cutpt18 = b;
        display('cutpt18 miss')
        display(i)
        end
        
        if Cutpt15(i) ~= min(Cutpt12(i),Cutpt15(i))
        b = [Cutpt15(1:i-1);0;Cutpt15(i:end)];
        Cutpt15 = b;
        display('cutpt15 miss')
        display(i)
        end
        
    end
end

Cutpt2dash = Cutpt2dashdash;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finished choosing the effective 2nd Cut point from Cutpt12, Cutpt15, Cutpt18 and
%Cuttpt2dash
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Section: 3   Receeding the chosen effective 2nd Cut point towards the collapse point
%for better estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Cutpt2 = zeros(length(Cutpt2dash),1);
%%%
%recede
L_array = Cutpt2dash-Cutpt1;
l_ini = max(L_array);
kappa = recede_fac_ini/power(l_ini,order);
for i=1:length(Cutpt2dash)
    recede_fac = 1-kappa*power(Cutpt2dash(i)-Cutpt1(i),order);
    Cutpt2(i) = floor((1-recede_fac)*Cutpt1(i)+recede_fac*Cutpt2dash(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finished receeding. Cutpt1 and Cutpt2 are now the two effecting cut
%points: One near the nucleation point and the other near the collapse
%point
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure
hold on
plot(t,RawSignal)
hold on
fplot(cutCurrentd,'-r')
hold on
fplot(cutCurrentr,'-g')
hold on
plot(t(Cutpt1),RawSignal(Cutpt1),'*r')
hold on
plot(t(Cutpt2),RawSignal(Cutpt2),'*k')
hold on
plot(t(Cutpt2),RawSignal(Cutpt2dash),'*g')
hold off
xlim([startpt endptt])
ylim([0 1000])
legend('Current','cutCurrentd','cutCurrentr','Cuptpt1','Cutpt2dash','Cutpt2', 'FontSize', 18)
ax = gca;
ax.FontSize = 18;
xlabel ('Time (\mus)', 'FontSize', 18)
ylabel ('Current (\muA)', 'FontSize', 18)


%%%cutting points%%%%%


% figure
% hold on
% plot(t,RawSignal)
% hold on
% plot(t(curvept1),RawSignal(curvept1),'*r')
% hold on
% plot(t(curvept2),RawSignal(curvept2),'*g')
% hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Section: 4   Start calculating the dip point from Cutpt1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
curvept1 = zeros(length(Cutpt1),1);
curvept2 = zeros(length(Cutpt2),1);

curvept1(1) = Cutpt1(1)-500;
for i = 1:length(Cutpt1)

    curvept2(i) = Cutpt1(i);
    curvept1(i+1) = Cutpt2(i);
end
dippt = zeros(length(curvept1)-1,1);
locationdippt = zeros(length(curvept1)-1,1);

% tolerance = 0.1;
factor = 0.005;

for i=1:length(curvept1)-1
% for i=1:1
    ptssavet = t(curvept1(i):curvept2(i));
%     if length(ptssavet)>50
        ptssaveRS = RawSignal(curvept1(i):curvept2(i));
%         f2 = fit(ptssavet,ptssaveRS,'fourier6');
%         ptsonf2 = f2(ptssavet);
        c = polyfit(ptssavet(length(ptssavet)-lengthfit(i):length(ptssavet)),ptssaveRS(length(ptssavet)-lengthfit(i):length(ptssavet)),1);
        
%         tolerance = factor * mean(abs(ptssavet(length(ptssavet)-lengthfit:length(ptssavet))-c(1)*ptssaveRS(length(ptssavet)-lengthfit:length(ptssavet))-c(2))/(sqrt(power(c(1),2)+1)));
%         tolerance = factor * abs(ptssavet(length(ptssavet))-c(1)*ptssaveRS(length(ptssavet))-c(2))/(sqrt(power(c(1),2)+1));
          tolerance = factor + abs((ptssaveRS(length(ptssavet))-c(2))/c(1)-ptssavet(length(ptssavet)));
%         d = zeros(length(ptssavet),1);
        d = tolerance - factor;
        j = length(ptssavet);
        marker = 1;
        
        disttt = zeros(length(ptssavet),1);
        while j >= 1 && marker == 1
%             d = abs(ptssavet(j)-c(1)*ptssaveRS(j)-c(2))/(sqrt(power(c(1),2)+1));
            d = abs((ptssaveRS(j)-c(2))/c(1)-ptssavet(j));
            disttt(j) = d;
            if d > tolerance
                dippt(i) = ptssavet(j-1);
                marker = 0;
            end
            j = j-1;
        end

           
        arr = find(t==dippt(i));
        locationdippt(i) = min(arr);
        
        
        if i<6 && i>2
        figure
        hold on
        plot(t,RawSignal)
        hold on
        plot(ptssavet,c(1)*ptssavet+c(2))
        hold on
        plot(ptssavet(length(ptssavet)),ptssaveRS(length(ptssavet)),'*b')
        hold on
        plot(dippt(i),RawSignal(locationdippt(i)),'*k')
        hold off
        hold on
        plot(t(Cutpt1),RawSignal(Cutpt1),'*r')
        hold on
        plot(t(Cutpt2),RawSignal(Cutpt2),'*g')
        hold off
        ylim([50 600])
        xlim([t(Cutpt2(i-1)) t(Cutpt1(i+1))])
        ylabel('Current (\muA)')
        xlabel('Time (\mus)')
        title('Dip point for 3rd bubble')
        end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of calculating the dip point from Cutpt1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% f1=figure
% % hold on
% % plot(t,RawSignal)
% % hold on
% [dippt_new,RS_loc_new]=dragpoints(t,RawSignal,dippt,RawSignal(locationdippt),min(t)-10,max(t)+10,min(RawSignal)-100,max(RawSignal)+300);
% % xlim([10 300])
% % hold off
% 
% % k = f1; %Some Figure
% %  while size(findobj(k))>0
% %     display('me'); %some action
% %     pause %some input
% %  end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Section: 5   Start calculating the rise point from Cutpt2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

risept = zeros(length(curvept1)-1,1);
locationrisept = zeros(length(curvept1)-1,1);

factor1 = 0.001;

for i=1:length(curvept1)-1
        ptssavet = t(curvept2(i):curvept1(i+1));
        ptssaveRS = RawSignal(curvept2(i):curvept1(i+1));

        lengthfit1 = averagingfac*length(ptssavet);
        c = polyfit(ptssavet(length(ptssavet)-lengthfit1:length(ptssavet)),ptssaveRS(length(ptssavet)-lengthfit1:length(ptssavet)),1);
        
        tolerance = factor + abs((ptssaveRS(length(ptssavet))-c(2))/c(1)-ptssavet(length(ptssavet)));

        d = tolerance - factor;
        j = length(ptssavet);
        marker = 1;
        
        disttt = zeros(length(ptssavet),1);
        while j >= 1 && marker == 1
            d = abs((ptssaveRS(j)-c(2))/c(1)-ptssavet(j));
            disttt(j) = d;
            if d > tolerance
                risept(i) = ptssavet(j+1);
                marker = 0;
            end
            j = j-1;
        end

    
        arr = find(t==risept(i));
        locationrisept(i) = min(arr);
        if i<4 && i>2
        figure
        hold on
        plot(t,RawSignal)
        hold on
        plot(ptssavet,c(1)*ptssavet+c(2))
        hold on
        plot(ptssavet(length(ptssavet)),ptssaveRS(length(ptssavet)),'*b')
        hold on
        plot(risept(i),RawSignal(locationrisept(i)),'*k')
        hold off
        hold on
        plot(t(Cutpt1),RawSignal(Cutpt1),'*r')
        hold on
        plot(t(Cutpt2),RawSignal(Cutpt2),'*g')
        hold off
        ylim([50 250])
        xlim([t(Cutpt2(i-1)) t(Cutpt1(i+1))])
        ylabel('Current (\muA)')
        xlabel('Time (\mus)')
        title('Rise point for 3rd bubble')
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%End of calculating the rise point from Cutpt2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



bubble_nos = length(curvept1)-1;
display(bubble_nos);

if length(dippt)~=length(risept)
    display("not ok")
end











waiting_time = zeros(length(curvept1)-1,1);
blockage_duration = zeros(length(curvept1)-1,1);
segment_nos = 10;
current_rise = zeros(length(curvept1)-1,segment_nos);
dipCurrent = zeros(length(curvept1)-1,1);
riseCurrent = zeros(length(curvept1)-1,1);
for i=1:length(curvept1)-1
    if i == 1
        waiting_time(i) = dippt(i)-startpt;
        dipCurrent(i) = mean(LowFilteredSignal(locationdippt(i)-nucpoint_averaging_width:locationdippt(i)+nucpoint_averaging_width));
        riseCurrent(i) = mean(LowFilteredSignal(locationrisept(i)-nucpoint_averaging_width:locationrisept(i)+nucpoint_averaging_width));

    else
        waiting_time(i) = dippt(i)-risept(i-1);
        dipCurrent(i) = mean(LowFilteredSignal(locationdippt(i)-nucpoint_averaging_width:locationdippt(i)+nucpoint_averaging_width));
        riseCurrent(i) = mean(LowFilteredSignal(locationrisept(i)-55:locationrisept(i)+nucpoint_averaging_width));
    end
    blockage_duration(i) = risept(i)-dippt(i);
end

locationdipptstart = locationdippt-nucpoint_averaging_width;
locationdipptend = locationdippt+nucpoint_averaging_width;





figure
hold on
plot(t,RawSignal)
hold on
plot(dippt,RawSignal(locationdippt),'*r')
hold on
plot(risept,RawSignal(locationrisept),'*g')
hold on
xlabel ('Time (\mus)')
ylabel ('Current (\muA)')
legend('Current','Nucleation point','Bubble collapse point')
% xlim([10 300])
hold off




figure
hold on
plot(t,LowFilteredSignal)
hold on
plot(dippt,LowFilteredSignal(locationdippt),'*r')
hold on

plot(t(locationdipptstart),LowFilteredSignal(locationdipptstart),'xr')
hold on
plot(t(locationdipptend),LowFilteredSignal(locationdipptend),'xr')
hold on
plot(risept,LowFilteredSignal(locationrisept),'*g')
hold on
xlabel ('Time (\mus)')
ylabel ('Current (\muA)')
legend('Current','Nucleation point','Bubble collapse point')
% xlim([10 300])
hold off




figure
histogram(blockage_duration(1:end), 40)
xlabel ('Blockage duration (\mus)')
ylabel ('Bubble number')
hold off

interestpts = vertcat(dippt,risept);
interestpts = sort(interestpts);
