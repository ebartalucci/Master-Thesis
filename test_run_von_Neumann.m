clear options

samples = sin(linspace(0,pi,181)).*linspace(0,pi,181);
samples = samples/sum(samples);

target = sin(linspace(0,pi,181)).^2;
target = target/sum(target);

options.trials = 10000;

tic, % starts timer
[distr,overlap,srate] = test_rejection_sampling(samples,target,options);
toc, % stops timer

fprintf(1,'Overlap = %4.3f\n',overlap);
fprintf(1,'Success rate = %8.5f\n',srate);

figure(1); clf; hold on;
plot(samples,'-','Color',[0.25,0.25,0.25]);
plot(target,'-','Color',[0,0.6,0]);
plot(distr,'-','Color',[0.75,0,0]);

M = max(target)/max(samples);
figure(2); clf; hold on;
plot(M*samples,'-','Color',[0.25,0.25,0.25]);
plot(target,'-','Color',[0,0.6,0]);
title('Scaling');

