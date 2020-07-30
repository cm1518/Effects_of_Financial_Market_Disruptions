shocks = [];
load basic_60bounds

r = randi([7500 15000],1, 500);
add_matrices2 = add_matrices(:,:,r);

%for i = 1:size(add_matrices2,3);
%    shocks(4*i+1-4:4*i,:) = add_matrices2(:,:,i);
%end

%shocks = shocks'max(;

%% Alternative version

for i = 1:size(add_matrices2,3);
    shock1(i,:) = add_matrices2(1,:,i);
    shock2(i,:) = add_matrices2(2,:,i);
    shock3(i,:) = add_matrices2(3,:,i);
    shock4(i,:) = add_matrices2(4,:,i);
end

shocks = [shock1; shock2; shock3; shock4];
shock3 = shock3(:,121:end)';
