%transition probability matrix

%blockage_duration_array


bins = 12;
class_width = 112;
class_start = 56;
trans_class = zeros(bins,bins);

%classification of bubbles
for counterpulse = 1:20
blockage_duration_pulse=blockage_duration_array(:,counterpulse);
blockage_duration_pulse(blockage_duration_pulse==0) = [];
class_dk1 = zeros(length(blockage_duration_pulse)-1,1);
class_prev_dk1 = zeros(length(blockage_duration_pulse)-1,1);

for i=2:length(blockage_duration_pulse)
    class_dk1(i-1) = ceil((1000*blockage_duration_pulse(i)-class_start)/class_width);
    class_prev_dk1(i-1) = ceil((1000*blockage_duration_pulse(i-1)-class_start)/class_width);
    trans_class(class_dk1(i-1),class_prev_dk1(i-1))=trans_class(class_dk1(i-1),class_prev_dk1(i-1))+1;
end

end

% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% class_dk2 = zeros(length(blokage_dk2),1);
% class_prev_dk2 = zeros(length(blokage_dk2)-1,1);
% for i=2:length(blokage_dk2)
%     class_dk2(i-1) = ceil((1000*blokage_dk2(i)-class_start)/class_width);
%     class_prev_dk2(i-1) = ceil((1000*blokage_dk2(i-1)-class_start)/class_width);
%     trans_class(class_dk2(i-1),class_prev_dk2(i-1))=trans_class(class_dk2(i-1),class_prev_dk2(i-1))+1;
% end
% class_dk3 = zeros(length(blokage_dk3),1);
% class_prev_dk3 = zeros(length(blokage_dk3)-1,1);
% for i=2:length(blokage_dk3)
%     class_dk3(i-1) = ceil((1000*blokage_dk3(i)-class_start)/class_width);
%     class_prev_dk3(i-1) = ceil((1000*blokage_dk3(i-1)-class_start)/class_width);
%     trans_class(class_dk3(i-1),class_prev_dk3(i-1))=trans_class(class_dk3(i-1),class_prev_dk3(i-1))+1;   
% end
% class_dk4 = zeros(length(blokage_dk4),1);
% class_prev_dk4 = zeros(length(blokage_dk4)-1,1);
% for i=2:length(blokage_dk4)
%     class_dk4(i-1) = ceil((1000*blokage_dk4(i)-class_start)/class_width);
%     class_prev_dk4(i-1) = ceil((1000*blokage_dk4(i-1)-class_start)/class_width);
%     trans_class(class_dk4(i-1),class_prev_dk4(i-1))=trans_class(class_dk4(i-1),class_prev_dk4(i-1))+1;    
% end
% 
% class_dk5 = zeros(length(blokage_dk5),1);
% class_prev_dk5 = zeros(length(blokage_dk5)-1,1);
% for i=2:length(blokage_dk5)
%     class_dk5(i-1) = ceil((1000*blokage_dk5(i)-class_start)/class_width);
%     class_prev_dk5(i-1) = ceil((1000*blokage_dk5(i-1)-class_start)/class_width);
%     trans_class(class_dk5(i-1),class_prev_dk5(i-1))=trans_class(class_dk5(i-1),class_prev_dk5(i-1))+1;    
% end
% 
% class_dk6 = zeros(length(blokage_dk6),1);
% class_prev_dk6 = zeros(length(blokage_dk6)-1,1);
% for i=2:length(blokage_dk6)
%     class_dk6(i-1) = ceil((1000*blokage_dk6(i)-class_start)/class_width);
%     class_prev_dk6(i-1) = ceil((1000*blokage_dk6(i-1)-class_start)/class_width);
%     trans_class(class_dk6(i-1),class_prev_dk6(i-1))=trans_class(class_dk6(i-1),class_prev_dk6(i-1))+1;    
% end
% 
% class_dk7 = zeros(length(blokage_dk7),1);
% class_prev_dk7 = zeros(length(blokage_dk7)-1,1);
% for i=2:length(blokage_dk7)
%     class_dk7(i-1) = ceil((1000*blokage_dk7(i)-class_start)/class_width);
%     class_prev_dk7(i-1) = ceil((1000*blokage_dk7(i-1)-class_start)/class_width);
%     trans_class(class_dk7(i-1),class_prev_dk7(i-1))=trans_class(class_dk7(i-1),class_prev_dk7(i-1))+1;    
% end
% 
% class_dk8 = zeros(length(blokage_dk8),1);
% class_prev_dk8 = zeros(length(blokage_dk8)-1,1);
% for i=2:length(blokage_dk8)
%     class_dk8(i-1) = ceil((1000*blokage_dk8(i)-class_start)/class_width);
%     class_prev_dk8(i-1) = ceil((1000*blokage_dk8(i-1)-class_start)/class_width);
%     trans_class(class_dk8(i-1),class_prev_dk8(i-1))=trans_class(class_dk8(i-1),class_prev_dk8(i-1))+1;    
% end
% 
% class_dk9 = zeros(length(blokage_dk9),1);
% class_prev_dk9 = zeros(length(blokage_dk9)-1,1);
% for i=2:length(blokage_dk9)
%     class_dk9(i-1) = ceil((1000*blokage_dk9(i)-class_start)/class_width);
%     class_prev_dk9(i-1) = ceil((1000*blokage_dk9(i-1)-class_start)/class_width);
%     trans_class(class_dk9(i-1),class_prev_dk9(i-1))=trans_class(class_dk9(i-1),class_prev_dk9(i-1))+1;    
% end
% 
% class_dk10 = zeros(length(blokage_dk10),1);
% class_prev_dk10 = zeros(length(blokage_dk10)-1,1);
% for i=2:length(blokage_dk10)
%     class_dk10(i-1) = ceil((1000*blokage_dk10(i)-class_start)/class_width);
%     class_prev_dk10(i-1) = ceil((1000*blokage_dk10(i-1)-class_start)/class_width);
%     trans_class(class_dk10(i-1),class_prev_dk10(i-1))=trans_class(class_dk10(i-1),class_prev_dk10(i-1))+1;    
% end

trans_prob = zeros(bins,bins);
for i=1:bins
    for j=1:bins
        trans_prob(i,j) = trans_class(i,j)/sum(trans_class(i,:));
    end
end
% figure
% for i=1:15
%     subplot(5,3,i)
%     bar(trans_prob(i,:))
%     title(['C',num2str(i),' #',num2str(sum(trans_class(i,:)))])
%     ylabel('P_t')
%     xticklabels({'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15'});
%     xtickangle(45)
% end
figure
heatmap(trans_prob);