x = [];

for i = 1:length(y)/4
   
    x(i) = sum(y(4*(i-1)+1:4*i))/4;
    
end

x = x';