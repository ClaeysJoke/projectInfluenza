data_week(:,1) = [41;42;43;44;45;46;47;48;49;50;51;0;1;2;3;4;5];
for j=2:5
for i=1:11   
    data_week(i,j) = ((data(2*i-1,j)*abs(data(2*i,1) - data_week(i,1))) + (data(2*i,j)*abs(data(2*i-1,1) - data_week(i,1))))/(data(2*i,1)-data(2*i-1,1));
end
for i=13:17   
    data_week(i,j) = ((data(2*i-1,j)*abs(data(2*i,1) - data_week(i,1))) + (data(2*i,j)*abs(data(2*i-1,1) - data_week(i,1))))/(data(2*i,1)-data(2*i-1,1));
end
data_week(12,j) = ((data(2*12-1,j)*abs(data(2*12,1) - data_week(12,1))) + (data(2*12,j)*abs(abs(52-data(2*12-1,1)) - data_week(12,1))))/(abs(data(2*12,1)+abs(data(2*12-1,1)-52)));
end

plot(data_week(:,2))
dlmwrite('data_week.txt', data_week, 'delimiter','\t','newline','pc')