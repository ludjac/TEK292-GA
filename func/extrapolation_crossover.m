function o = extrapolation_crossover(p1, p2, fitfunc, interval)
c = 0.5;
o1 = [c*p1(3) + c*p2(3), 3*c*p1(3) - c*p2(3), -c*p1(3) + 3*c*p2(3)];
o2 = [c*p1(4) + c*p2(4), 3*c*p1(4) - c*p2(4), -c*p1(4) + 3*c*p2(4)];
index1 = find(abs(o1)<interval(2));
index2 = find(abs(o2)<interval(2));
o1 = o1(index1);
o2 = o2(index2);

fitness = [];

for i=1:1:length(o1)
    for k=1:1:length(o2)
       fitness = [fitness; i k fitfunc(o1(i), o2(k))];
    end
end
% Minimize (2) / Maximize (2)
fitness = sortrows(fitness, [2]);
o = [o1(fitness(1, 1)) o2(fitness(1,2))];
