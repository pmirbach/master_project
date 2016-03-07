function [h] = test2()

h(1) = figure;
subplot(1,2,1)
plot(1:10)
legend('blub')
subplot(1,2,2)
plot(10:-1:1)

h(2) = figure;
subplot(1,2,1)
plot(1:10)
subplot(1,2,2)
plot(10:-1:1)