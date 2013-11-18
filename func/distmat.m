function p2 = distmat(p, pop)

[len dim] = size(pop);

D = NaN(len-1, 2);
for i=2:1:len
    D(i-1, :) = [pop(i,1), sqrt(sum(((p(3:4)-pop(i,3:4)).^2)))];
end
D_sorted = sortrows(D, [2]);
thresh = ceil(len*0.25);
index = randi([1, thresh], 1);
ind = D_sorted(index, 1); 
p2 = pop(pop(:, 1)==ind, :);