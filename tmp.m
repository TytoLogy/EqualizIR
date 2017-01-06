load s

figure(1)

plot(1:512, s(1:512), 1:512, fliplr(s(513:end)))

sum(s(1:512) - fliplr(s(513:end)))

