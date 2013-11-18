function ind = GAwheel(pop)

pop = sortrows(pop, [2]);
s = sum(pop(:,2));

prob = [pop(:,2)./s];

c = cumsum(prob);

i = 0;
r = rand(1);
temp = find((c < r));
if isempty(temp)
    i = 1;
else
    i = temp(end)+1;
end

ind = pop(i, :);
