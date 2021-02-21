% this script generates distance distributions from an ensemble of hnRNP A1
% conformers by calling rotamer library simulation in MMMx

pairs = zeros(19,2); % initialize array that specifies the site pairs 
pairs(1,:) = [52,231]; % please check against hnRNPA1_flex_restraints.mcx
pairs(2,:) = [144,231];
pairs(3,:) = [182,231];
pairs(4,:) = [52,271];
pairs(5,:) = [144,271];
pairs(6,:) = [182,271];
pairs(7,:) = [52,316];
pairs(8,:) = [144,316];
pairs(9,:) = [182,316];
pairs(10,:) = [32,231];
pairs(11,:) = [52,197];
pairs(12,:) = [182,190];
pairs(13,:) = [182,197];
pairs(14,:) = [182,223];
pairs(15,:) = [182,252];
pairs(16,:) = [182,297];
pairs(17,:) = [231,271];
pairs(18,:) = [231,316];
pairs(19,:) = [271,316];

% Determine number m of restraints
[m,~] = size(pairs);

% find all conformers in the current directory
all_files = dir('hnRNPA1_NMR_unrestrained_*.pdb');

% read the first file in order to initialize distance axis and fit its
% length
entity = get_pdb(all_files(1).name); % this is the PDB file reader of MMMx
site1 = sprintf('(A)%i',pairs(1,1));
site2 = sprintf('(A)%i',pairs(1,2));
% the following is how MMMx generates a distance distribution
[r_axis,distribution,entity,exceptions] = distance_distribution(entity,site1,'mtsl',site2,'mtsl');

% plot it to check that everything is right
figure(1); clf;
plot(r_axis,distribution);

all_distributions = zeros(m,length(r_axis));

tic,
for k = 1:length(all_files) % loop over all conformers
    entity = get_pdb(all_files(k).name);
    % Keep the users happy by telling them that something happens 
    fprintf(1,'Processing %s\n',all_files(k).name);
    for kr = 1:m % loop over all site pairs
        site1 = sprintf('(A)%i',pairs(kr,1)); % first site address
        site2 = sprintf('(A)%i',pairs(kr,2)); % second site address
        fprintf(1,'   %s-%s\n',site1,site2);
        % get distance distribution, distance axis will be the same
        % get entity (protein structure) back from the MMMx routine, this
        % saves computation effort when later using one of the labeling
        % sites again
        [~,distribution,entity] = distance_distribution(entity,site1,'mtsl',site2,'mtsl');
        % sometimes in-silico labeling fails, we have to guard against that
        if ~isempty(distribution) % only if labeling was successful
            % add current distribution to ensemble distribution
            all_distributions(kr,:) = all_distributions(kr,:) + distribution;
        end
    end
end
toc,

% open a file for writing that reports mean values and standard deviations
fid = fopen('sim_restraints.dat','wt');
for kr = 1:m % loop over all restraints (site pairs)
    % plot ensemble distance distribution
    figure(kr); clf;
    distribution = all_distributions(kr,:);
    distribution = distribution/sum(distribution); % it should be normalized
    plot(r_axis,distribution);
    % make a nice axis scaling, 0... 100 Angstroem, vertical axis not
    % completely tight
    axis([10,100,-0.1*max(distribution),1.1*max(distribution)]);
    % convert from Angstroem to nanometers and add artificial lower and
    % upper bounds
    % as we have generated row vectors, we need to convert to column
    % vectors by the transpose (.')
    data = [r_axis.'/10 distribution.' 0.9*distribution.' 1.1*distribution.'];
    % make the file name for the distance distribution file of this site
    % pair
    fname = sprintf('sim-%i-%i-distr.dat',pairs(kr,1),pairs(kr,2));
    % use the file name as figure title
    title(fname);
    save(fname,'data','-ascii');
    % compute the mean value
    rmean = sum(r_axis.*distribution);
    % compute standard deviation
    stdr = sqrt(sum(distribution.*(r_axis-rmean).^2));
    % print site addresses, mean value, and standard deviation to restraint
    % file
    fprintf(fid,'  (A)%i   (A)%i  %4.2f  %4.2f\n',pairs(kr,1),pairs(kr,2),rmean,stdr);
end
fclose(fid);