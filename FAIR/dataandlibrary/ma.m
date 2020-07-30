x = [];

for i = 1:length(y)/12
   
    x(i) = sum(y(12*(i-1)+1:12*i))/12;
    
end

x = x';