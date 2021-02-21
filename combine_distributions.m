data = load('DEER_A1v33R1_K52C_S231C_DeerLab_alpha_5_distr.dat');
r_axis = 10*data(:,1); % convert to Angstroem
% distr = data(:,2);
% distr = distr/sum(distr); % normalize

% or take the lower bound
distr = data(:,3);
distr = distr/sum(distr); % normalize

% ### This needs to be corrected ###
% data2 = load('sim-52-231-distr.dat');
data2 = load('sim_190_231_distr.dat'); 
r_axis2 = 10*data2(:,1);
distr2 = data2(:,2);

% we interpolate to the distance axis of the experimental data
distr2b = interp1(r_axis2,distr2,r_axis,'pchip',0);
distr2b = distr2b/sum(distr2b); % normalize

% combine the data vectors
% first column is distance axis, second column is g(r), third column is
% f(r)
data = [r_axis distr2b distr];

save('input-vNrs-52-321.dat','data','-ascii');

% for reading, use
data = load('input-vNrs-52-321.dat');
r_axis = data(:,1);
samples = data(:,2);
target = data(:,3);

% check that it worked

figure(1); clf; hold on
plot(r_axis,distr,'k');
plot(r_axis,distr2b,'r');
xlabel('distance (nm)');
ylabel('Probability density (a.u.)');

