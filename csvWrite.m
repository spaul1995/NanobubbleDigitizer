function M = csvWrite(nos,pulsewidth,t1,Iq,Iq2)

starttimept1 = 2*pulsewidth*(nos-1)+1;
endtimept1 = 2*pulsewidth*(nos-1)+(pulsewidth-0.5);
startindex = length(t1(t1<starttimept1));
endindex = length(t1(t1<endtimept1));

M(1:endindex-startindex+1,1) = t1(startindex:endindex);
M(1:endindex-startindex+1,2) = Iq(startindex:endindex);
M(1:endindex-startindex+1,3) = Iq2(startindex:endindex);
% writematrix('M2.csv',M2)
filename = sprintf('M%i.csv',fix(nos));
dlmwrite(filename, M, 'precision', '%9i')
end

