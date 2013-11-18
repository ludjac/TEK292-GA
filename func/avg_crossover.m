function o = avg_crossover(p1, p2, ~, ~)
avg = @(p1, p2) (p1 + p2)/2;
o = [avg(p1(3), p2(3)) avg(p1(4), p2(4))];