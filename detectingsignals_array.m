% NanobubbleDigitizer data analysis package licence details:
% 
% Copyright Soumyadeep Paul (soumyadeep.paul@thml.t.u-tokyo.ac.jp)
% 
% This file is part of NanobubbleDigitizer data analysis package.
% 
% NanobubbleDigitizer is a MATLAB package designed to detect discrete 
% bubble events from transient current traces. The bubble features
% extracted from the current data using this code can be used to analyze
% nanoscopic bubble dynamics confined within solid-state nanopores. This
% package is released under the GNU GPL.
% 
% The copyright of arb is held by Soumyadeep Paul.
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
% For full details of arb's licence see the licence file in the main 
% directory.
clear all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pulses = 20;
endingarray = [246250,246250,246250,246250,243250,246250,246250,246250,246250,245250,245000,246250,246250,246250,246250,244250,244250,246250,246250,235550];

startptarray = [0,200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200,3400,3600,3800];
endpttarray = startptarray+100;

base1array = [220,220,220,220,220,220,220,220,220,220,220,220,220,220,220,220,220,220,220,220];
base2array = [310,290,290,290,290,290,290,290,290,290,290,290,290,290,290,290,290,290,290,290];

recede_fac_iniarray = [0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03];
orderarray = [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5];

averagingfacarray=[20/1250,20/1250,20/1250,20/1250,20/1250,20/1250,20/1250,20/1250,20/1250,20/1250,20/1250,20/1250,20/1250,20/1250,20/1250,20/1250,20/1250,20/1250,20/1250,20/1250];

nucpoint_averaging_widtharray = fix((28.* ones(20,1))');

waiting_time_array = zeros(100,pulses);
blockage_duration_array = zeros(100,pulses);
dipCurrent_array = zeros(100,pulses);
riseCurrent_array = zeros(100,pulses);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for counterpulse = 1:20
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startpt = startptarray(counterpulse);
endptt = endpttarray(counterpulse);

% base1 = base1array(counterpulse);
% base2 = base2array(counterpulse);

base1 = 180;
base2 = 350;

base12 = 200;
base15 = 250;
base18 = 300;

recede_fac_ini = recede_fac_iniarray(counterpulse);
order = orderarray(counterpulse);

averagingfac=averagingfacarray(counterpulse);

nucpoint_averaging_width = nucpoint_averaging_widtharray(counterpulse);

ending = endingarray(counterpulse);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 



filename = sprintf('M%i.csv',fix(counterpulse));
display(filename)
start = 0;
%ending = 137600;

%ending = 621150;
X1 = csvread(filename,start,0,[start,0,ending,2]);

global V;
global pr;
global alpha_sigma;
V = 10;
alpha_sigma = 2.52;
global room_temperature;
room_temperature = 300.15;
pr = 500;
%control parameters

% T = 0.0008*factor;
T1 = 0.0004;
RawSignal = X1(:,2);
LowFilteredSignal = X1(:,3);
% RawSignal1 = K_1;
% L = size(RawSignal);
% L1 = size(RawSignal1);
t = X1(:,1);
% t1 = (0:L1-1)*T1;
StartCoeff = 2;
EndCoeff = 0;
FilterCoeff = 0.999;
% cutCurrent = 120;
K = 10e5;
s1 = 950;
s2 = s1-450;
cutCurrentd = @(t) s1-base1-s2*(1-exp(-K*1e-6*(t-startpt)));

cutCurrent12 = @(t) s1-base12-s2*(1-exp(-K*1e-6*(t-startpt)));
cutCurrent15 = @(t) s1-base15-s2*(1-exp(-K*1e-6*(t-startpt)));
cutCurrent18 = @(t) s1-base18-s2*(1-exp(-K*1e-6*(t-startpt)));

cutCurrentr = @(t) s1-base2-s2*(1-exp(-K*1e-6*(t-startpt)));
% figure
% fplot(cutCurrentd)
% xlim([3 50])
sigma = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% endpt = startpt+5;
% startindex = 1;
% endindex = length(t(t<endpt));
% meanAvg = mean(RawSignal(startindex:endindex));
% sigma = 5;
% Icut = RawSignal-cutCurrent;
% AbsIcut = abs(Icut);
% AbsIcut(AbsIcut>sigma) = cutCurrent;
% [pks,locs] = findpeaks(-AbsIcut);
% 
% 
% 
% bubble_number = 1;
% % saveindex1 = zeros(length(locs),1);
% % for i = 1:length(RawSignal)
% %     if RawSignal(i)<100+sigma && RawSignal(i)>100-sigma
% %         saveindex1(bubble_number) = i;
% %         bubble_number=bubble_number+1;
% %     end
% % end
% Cutpt1 = zeros(floor(length(locs)/2),1);
% Cutpt2 = zeros(floor(length(locs)/2),1);
% count = 1;
% switch1=1;
% marker = 0;
% 
% 
% for i = 1:length(RawSignal)
%     if (RawSignal(i)>cutCurrent+sigma)
%         marker = 0;
%     end
% if (RawSignal(i)<cutCurrent+sigma) && marker ==0
%     if ismember(i,locs)
%         if switch1==1
%         Cutpt1(bubble_number) = i;
%         switch1 = 2;
%         end
%     end
%     
% if switch1==2
%     bubble_number = bubble_number+1;
%     switch1 = 1;
%     marker = 1;
% end
% end
% end
% 
% Cutpt1(Cutpt1==0)=[];
% 
% 
% switch1=1;
% marker = 0;
% locs1 = setdiff(locs,Cutpt1);
% 
% for i =1:length(Cutpt1)-1
%     Cutpt2(i) = locs(length(locs(locs<Cutpt1(i+1))));
% end
% Cutpt2(i+1) = locs(length(locs));
% Cutpt2(Cutpt2==0)=[];



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
spacer = 20;
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
   
   
  if (RawSignal(i)<cutCurrentr(t(i))+sigma/2) && (RawSignal(i)>cutCurrentr(t(i))-sigma/2) && sign(RawSignal(i+1)-RawSignal(i)) > 0 && switch2 == 0
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

len_1 = sprintf('lenght of Cutpt1 = %i',fix(length(Cutpt1)));
len_12 = sprintf('lenght of Cutpt12 = %i',fix(length(Cutpt12)));
len_15 = sprintf('lenght of Cutpt15 = %i',fix(length(Cutpt15)));
len_18 = sprintf('lenght of Cutpt18 = %i',fix(length(Cutpt18)));
len_2 = sprintf('lenght of Cutpt2dash = %i',fix(length(Cutpt2dash)));
display(len_1)
display(len_12)
display(len_15)
display(len_18)
display(len_2)


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

% 
% 
% figure
% hold on
% plot(t,RawSignal)
% hold on
% fplot(cutCurrentd,'-r')
% hold on
% fplot(cutCurrentr,'-g')
% hold on
% plot(t(Cutpt1),RawSignal(Cutpt1),'*r')
% hold on
% plot(t(Cutpt2),RawSignal(Cutpt2),'*g')
% hold off
% xlim([startpt endptt])
% ylim([0 1000])

curvept1 = zeros(length(Cutpt1),1);
curvept2 = zeros(length(Cutpt2),1);

curvept1(1) = Cutpt1(1)-500;
for i = 1:length(Cutpt1)

    curvept2(i) = Cutpt1(i);
    curvept1(i+1) = Cutpt2(i);
end

%%%cutting points%%%%%


% figure
% hold on
% plot(t,RawSignal)
% hold on
% plot(t(curvept1),RawSignal(curvept1),'*r')
% hold on
% plot(t(curvept2),RawSignal(curvept2),'*g')
% hold off


dippt = zeros(length(curvept1)-1,1);
locationdippt = zeros(length(curvept1)-1,1);

% tolerance = 0.1;
factor = 0.01;

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

        
        


        
%         error = abs(ptsonf2-ptssaveRS);
%         
%         diagonaldist = sqrt(power(ptssavet,2)+power(ptsonf2,2));
%         [M1,I1] = max(diagonaldist);
%         dippt(i) = ptssavet(I1);
        
        
%         [d1,d2] = differentiate(f2,ptssavet);
%         [pks1,locs1] = findpeaks(-abs(d1));
%         dippt(i) = ptssavet(locs1(end));

%     else
%         [M,I] = max(RawSignal(curvept1(i):curvept2(i)));
%         dippt(i) = ptssavet(I);
%     end
        
        arr = find(t==dippt(i));
        locationdippt(i) = min(arr);
        
        
%         if i<4 && i>2
%         figure
%         hold on
%         plot(t,RawSignal)
%         hold on
%         plot(ptssavet,c(1)*ptssavet+c(2))
%         hold on
%         plot(ptssavet(length(ptssavet)),ptssaveRS(length(ptssavet)),'*b')
%         hold on
%         plot(dippt(i),RawSignal(locationdippt(i)),'*k')
%         hold off
%         hold on
%         plot(t(Cutpt1),RawSignal(Cutpt1),'*r')
%         hold on
%         plot(t(Cutpt2),RawSignal(Cutpt2),'*g')
%         hold off
%         ylim([50 250])
%         xlim([t(Cutpt2(i-1)) t(Cutpt1(i+1))])
%         ylabel('Current (\muA)')
%         xlabel('Time (\mus)')
%         title('Dip point for 3rd bubble')
%         end
end








risept = zeros(length(curvept1)-1,1);
locationrisept = zeros(length(curvept1)-1,1);

factor1 = 0.001;

for i=1:length(curvept1)-1
        ptssavet = t(curvept2(i):curvept1(i+1));
        ptssaveRS = RawSignal(curvept2(i):curvept1(i+1));
%         [posloc_min_val,posloc_min_index] = findpeaks(-ptssaveRS);
%         posloc_min = ptssavet(posloc_min_index);
% %         lengthfit1 = 0.7*(length(ptssavet)-posloc_min_index(end));
%         
%         laplacian = abs(gradient(gradient(ptssaveRS)));
%         [posloc_min_val1,posloc_min_index1] = findpeaks(-laplacian);
%         if isempty(posloc_min_val1) == 1
%            lengthfit1 = 0.7*(length(ptssavet)-posloc_min_index(end)) ;
%         else
%           posloc_min1 = ptssavet(posloc_min_index1);
%         lengthfit1 = (length(ptssavet)-max(posloc_min_index1(end), posloc_min_index(end)));  
%         end
        
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

        
%         
%         figure
%         hold on
%         plot(ptssavet,ptssaveRS)
%         hold on
%         plot(posloc_min, -posloc_min_val,'*g')
%         hold off

        
%         if i<3
%         figure
%         hold on
%         plot(ptssavet,ptssaveRS)
%         hold on
%         plot(ptssavet,c(1)*ptssavet+c(2))
%         hold on
%         plot(ptssavet(length(ptssavet)),ptssaveRS(length(ptssavet)),'*g')
%         hold on
%         plot(ptssavet(j+1),ptssaveRS(j+1),'*b')
%         hold off
%         ylim([0 200])
%         xlim([startpt endptt])
%         ylabel('Current (\muA)')
%         xlabel('Time (\mus)')
%         end

        arr = find(t==risept(i));
        locationrisept(i) = min(arr);
%         if i<4 && i>2
%         figure
%         hold on
%         plot(t,RawSignal)
%         hold on
%         plot(ptssavet,c(1)*ptssavet+c(2))
%         hold on
%         plot(ptssavet(length(ptssavet)),ptssaveRS(length(ptssavet)),'*b')
%         hold on
%         plot(risept(i),RawSignal(locationrisept(i)),'*k')
%         hold off
%         hold on
%         plot(t(Cutpt1),RawSignal(Cutpt1),'*r')
%         hold on
%         plot(t(Cutpt2),RawSignal(Cutpt2),'*g')
%         hold off
%         ylim([50 250])
%         xlim([t(Cutpt2(i-1)) t(Cutpt1(i+1))])
%         ylabel('Current (\muA)')
%         xlabel('Time (\mus)')
%         title('Rise point for 3rd bubble')
%         end
end


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

waiting_time_array(1:length(curvept1)-1,counterpulse) = waiting_time;
blockage_duration_array(1:length(curvept1)-1,counterpulse) = blockage_duration;
dipCurrent_array(1:length(curvept1)-1,counterpulse) = dipCurrent;
riseCurrent_array(1:length(curvept1)-1,counterpulse) = riseCurrent;




end
waiting_time_1Darray = waiting_time_array(:,1);
blockage_duration_1Darray = blockage_duration_array(:,1);
dipCurrent_1Darray = dipCurrent_array(:,1);
riseCurrent_1Darray = riseCurrent_array(:,1);

for i=2:20
waiting_time_1Darray = vertcat(waiting_time_1Darray,waiting_time_array(:,i));
blockage_duration_1Darray = vertcat(blockage_duration_1Darray,blockage_duration_array(:,i));
dipCurrent_1Darray = vertcat(dipCurrent_1Darray,dipCurrent_array(:,i));
riseCurrent_1Darray = vertcat(riseCurrent_1Darray,riseCurrent_array(:,i));
end
waiting_time_1Darray(waiting_time_1Darray==0) = [];
blockage_duration_1Darray(blockage_duration_1Darray==0) = [];
dipCurrent_1Darray(dipCurrent_1Darray==0) = [];
riseCurrent_1Darray(riseCurrent_1Darray==0) = [];

M(:,1) = 1000*waiting_time_1Darray;
M(:,2) = 1000*blockage_duration_1Darray;
M(:,3) = dipCurrent_1Darray;
M(:,4) = riseCurrent_1Darray;
M(:,5) = 9.0;
dlmwrite('9V_1micron_current.csv', M, 'precision', '%9i')


figure
g = histogram(1000*blockage_duration_1Darray);
g.BinWidth = 28;
set(g,'facecolor',[0.3010, 0.7450, 0.9330])
set(g,'edgecolor',[1, 1, 1])
set(gca,'FontSize',18)
% yt = get(gca, 'YTick');
% set(gca, 'YTick', yt, 'YTickLabel')
xlim([0 1500])
xlabel ('Blockage duration (ns)','FontSize', 20)
ylabel ('Frequency','FontSize', 20)
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
% for i =1:length(dippt)
% plot(dippt(i),RawSignal(locationdippt(i)),'*r')
% % txt = sprintf('Nuc @ %.2f \mus \rightarrow', dippt(i));
% % text(dippt(i),RawSignal(locationdippt(i)),txt,'HorizontalAlignment','right')
% end
% hold on
% for i =1:length(risept)
% plot(risept(i),RawSignal(locationrisept(i)),'*g')
% % txt = sprintf('Col @ %.2f \mus \rightarrow', risept(i));
% % text(risept(i),RawSignal(locationrisept(i)),txt,'HorizontalAlignment','right')
% end
hold on
plot(risept,LowFilteredSignal(locationrisept),'*g')
hold on
xlabel ('Time (\mus)')
ylabel ('Current (\muA)')
legend('Current','Nucleation point','Bubble collapse point')
% xlim([10 300])
hold off



% figure
% hold on
% plot(t,RawSignal)
% hold on
% plot(dippt,RawSignal(locationdippt),'*r')
% % for i =1:length(dippt)
% % plot(dippt(i),RawSignal(locationdippt(i)),'*r')
% % % txt = sprintf('Nuc @ %.2f \mus \rightarrow', dippt(i));
% % % text(dippt(i),RawSignal(locationdippt(i)),txt,'HorizontalAlignment','right')
% % end
% % hold on
% % for i =1:length(risept)
% % plot(risept(i),RawSignal(locationrisept(i)),'*g')
% % % txt = sprintf('Col @ %.2f \mus \rightarrow', risept(i));
% % % text(risept(i),RawSignal(locationrisept(i)),txt,'HorizontalAlignment','right')
% % end
% hold on
% plot(risept,RawSignal(locationrisept),'*g')
% hold on
% waiting_time = zeros(length(curvept1)-1,1);
% blockage_duration = zeros(length(curvept1)-1,1);
% segment_nos = 10;
% current_rise = zeros(length(curvept1)-1,segment_nos);
% dipCurrent = zeros(length(curvept1)-1,1);
% riseCurrent = zeros(length(curvept1)-1,1);
% for i=1:length(curvept1)-1
%     if i == 1
%         waiting_time(i) = dippt(i)-startpt;
%         dipCurrent(i) = RawSignal(locationdippt(i));
%         riseCurrent(i) = RawSignal(locationrisept(i));
%         X = RawSignal(1:locationdippt(i));
%         TT = t(1:locationdippt(i));
%         r = diff(fix(linspace(0, length(TT), segment_nos+1)));
%         display(length(TT))
%         current_subarray = mat2cell(X, r, 1);
%         time_subarray = mat2cell(TT, r, 1);
%         for j=1:segment_nos
%             p = polyfit(time_subarray{j,1},current_subarray{j,1},1);
%             f1 = polyval(p,time_subarray{j,1});
%             plot(time_subarray{j,1}(1),f1(1),'o')
%             hold on
%             plot(time_subarray{j,1}(end),f1(end),'o')
%             hold on
%             plot(time_subarray{j,1},f1,'r--')
%             hold on
%             current_rise(i,j) = p(1);
%         end
%        
%     else
%         waiting_time(i) = dippt(i)-risept(i-1);
%         dipCurrent(i) = RawSignal(locationdippt(i));
%         riseCurrent(i) = RawSignal(locationrisept(i));
%         X = RawSignal(locationrisept(i-1):locationdippt(i));
%         TT = t(locationrisept(i-1):locationdippt(i));
%         r = diff(fix(linspace(0, length(TT), segment_nos+1)));
%         display(length(TT))
%         current_subarray = mat2cell(X, r, 1);
%         time_subarray = mat2cell(TT, r, 1);
%         for j=1:segment_nos
%             p = polyfit(time_subarray{j,1},current_subarray{j,1},1);
%             f1 = polyval(p,time_subarray{j,1});
%             plot(time_subarray{j,1}(1),f1(1),'o')
%             hold on
%             plot(time_subarray{j,1}(end),f1(end),'o')
%             hold on
%             plot(time_subarray{j,1},f1,'r--')
%             hold on
%             current_rise(i,j) = p(1);
%         end
%     end
%     blockage_duration(i) = risept(i)-dippt(i);
% end
% 
% 
% 
% xlabel ('Time (\mus)')
% ylabel ('Current (\muA)')
% legend('Current','Nucleation point','Bubble collapse point')
% % xlim([10 300])
% hold off
% 
% 
% 
% 
% 
% figure
% hold on
% plot(t,RawSignal)
% hold on
% plot(dippt,RawSignal(locationdippt),'*r')
% % for i =1:length(dippt)
% % plot(dippt(i),RawSignal(locationdippt(i)),'*r')
% % % txt = sprintf('Nuc @ %.2f \mus \rightarrow', dippt(i));
% % % text(dippt(i),RawSignal(locationdippt(i)),txt,'HorizontalAlignment','right')
% % end
% % hold on
% % for i =1:length(risept)
% % plot(risept(i),RawSignal(locationrisept(i)),'*g')
% % % txt = sprintf('Col @ %.2f \mus \rightarrow', risept(i));
% % % text(risept(i),RawSignal(locationrisept(i)),txt,'HorizontalAlignment','right')
% % end
% hold on
% plot(risept,RawSignal(locationrisept),'*g')
% hold on
% 
% 
% 
% segment_nos = 15;
% current_rise2 = zeros(length(curvept1)-1,segment_nos);
% for i=1:length(curvept1)-1
%     if i == 1
%         TT = t(1:locationdippt(i));
%         r = diff(fix(linspace(0, length(TT), segment_nos+1)));
%         r = cumsum(r);
%         r2 = [1,r];
%         X = zeros(length(r)+1,1);
%         TT = zeros(length(r)+1,1);
%         for k=1:length(r2)
%             X(k) = RawSignal(r2(k));
%             TT(k) = t(r2(k));
%         end
%         for j=1:segment_nos
%             plot(TT(j),X(j),'o')
%             hold on
%             plot(linspace(TT(j),TT(j+1)),linspace(X(j),X(j+1)),'r--')
%             hold on
%             current_rise2(i,j) = (X(j+1)-X(j))/(TT(j+1)-TT(j));
%         end
%        plot(TT(j+1),X(j+1),'o')
%     else
%         TT = t(locationrisept(i-1):locationdippt(i));
%         r = diff(fix(linspace(0, length(TT), segment_nos+1)));
%         r = locationrisept(i-1)+cumsum(r);
%         r2 = [locationrisept(i-1),r];
%         X = zeros(length(r)+1,1);
%         TT = zeros(length(r)+1,1);
%         for k=1:length(r2)
%             X(k) = RawSignal(r2(k));
%             TT(k) = t(r2(k));
%         end
%         for j=1:segment_nos
%             plot(TT(j),X(j),'o')
%             hold on
%             plot(linspace(TT(j),TT(j+1)),linspace(X(j),X(j+1)),'r--')
%             hold on
%             current_rise2(i,j) = (X(j+1)-X(j))/(TT(j+1)-TT(j));
%         end
%         plot(TT(j+1),X(j+1),'o')
%     end
% end
% 
% 
% 
% xlabel ('Time (\mus)')
% ylabel ('Current (\muA)')
% legend('Current','Nucleation point','Bubble collapse point')
% % xlim([10 300])
% hold off



figure
histogram(blockage_duration(1:end), 40)
xlabel ('Blockage duration (\mus)')
ylabel ('Bubble number')
hold off
% 
% figure
% histogram(waiting_time(2:end), 10)
% xlabel ('Waiting time (\mus)')
% ylabel ('Bubble number')
% hold off
% 
% figure
% hold on
% yyaxis left
% plot(waiting_time(2:end),'--*r')
% ylabel('Waiting time (\mus)')
% hold on
% yyaxis right
% plot(blockage_duration(2:end),'--*g')
% ylabel('Blockage duration (\mus)')
% xlabel('Bubble number')
% legend('Waiting time','Blockage duration')
% hold off
% 
% figure
% bar(1000*blockage_duration,'facecolor','r')
% ylabel('Blockage duration (ns)')
% xlabel('Bubble number')
% 
% figure
% bar(waiting_time(1:end),'facecolor','r')
% ylabel('Waiting time (\mus)')
% xlabel('Bubble number')
% ylim([0 10])

interestpts = vertcat(dippt,risept);
interestpts = sort(interestpts);
% tic
% transientJH_wall;
% toc

% 
% figure
% hold on
% yyaxis left
% plot(ptssavet,ptssaveRS)
% hold on
% yyaxis left
% plot(f2,ptssavet,ptssaveRS)
% hold on
% yyaxis right
% plot(ptssavet,d1)
% hold on
% yyaxis right
% plot(ptssavet(locs1(end)),d1(locs1(end)),'*g')
% hold off




