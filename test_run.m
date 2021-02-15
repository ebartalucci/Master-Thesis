clear options

samples = sin(linspace(0,pi,181)).^2;

options.trials = 10000;

tic, % starts timer
[distr,overlap] = test_rejection_sampling(samples,[],options);
toc, % stops timer

fprintf(1,'Overlap = %4.3f\n',overlap);

figure(1); clf; hold on;
plot(samples/sum(samples),'-','Color',[0.25,0.25,0.25]);
plot(distr,'-','Color',[0.75,0,0]);
