tic;
% *1. Data preprocessing steps*
%% 
% _*(1.1)* Load the fMRI data for both control and patients_

%read in the bipolar all at once
bipolar_timeseries = {};
files = dir(fullfile('/home/tonnar/Fall2021/Class/BME/FinalProject/Bipolar_Task', "*.csv"));
for i = 1:length(files)
    bipolar_ts = readmatrix(files(i).name);
    bipolar_ts(:,1) = [];
    bipolar_ts(1,:) = [];
    bipolar_timeseries{i} = bipolar_ts';
end
clear bipolar_ts

%read in the control all at once
control_timeseries = {};
files = dir(fullfile('/home/tonnar/Fall2021/Class/BME/FinalProject/Control_Task', "*.csv"));
for i = 1:length(files)
    control_ts = readmatrix(files(i).name);
    control_ts(:,1) = [];
    control_ts(1,:) = [];
    control_timeseries{i} = control_ts';
end
clear control_ts
clear files
%% 
% _*(1.2)* Normalize the data and calculate correlation matricies_

%normalize the data
for i = 1:length(bipolar_timeseries)
    bipolar_timeseries{i} = bipolar_timeseries{i} - mean(bipolar_timeseries{i}, 2,"omitnan");
    bipolar_timeseries{i} = bipolar_timeseries{i} ./ vecnorm(bipolar_timeseries{i}, 2, 2);
end

for i = 1:length(control_timeseries)
    control_timeseries{i} = control_timeseries{i} - mean(control_timeseries{i}, 2,"omitnan");
    control_timeseries{i} = control_timeseries{i} ./ vecnorm(control_timeseries{i}, 2, 2);
end

%check if it was normalized correctly
means = [];
for i = 1:length(bipolar_timeseries)
    means = [means, max(mean(bipolar_timeseries{i}, 2,"omitnan"))];
    means = [means, max(mean(control_timeseries{i}, 2,"omitnan"))];
end

norms = [];
for i = 1:length(bipolar_timeseries)
    norms = [norms, vecnorm(control_timeseries{i}, 2, 2)];
    norms = [norms, vecnorm(control_timeseries{i}, 2, 2)];
end

disp("The max mean: " + round(max(means), 10) + " The min mean: " + round(min(means), 10));
disp("The max norm: " + max(max(norms)) + " The min norm: " + min(min(norms)))
%calculate the correlations
bipolar_corr = {};
for i = 1:length(bipolar_timeseries)
    bipolar_corr{i} = corr(bipolar_timeseries{i}');
end

control_corr = {};
for i = 1:length(control_timeseries)
    control_corr{i} = corr(control_timeseries{i}');
end
%%

%lets visualize 
for i = 1:length(bipolar_timeseries)
    vis_timeseries(bipolar_timeseries{i})
end

for i = 1:length(control_timeseries)
    vis_timeseries(control_timeseries{i})
end
% *2. Network Construction*
%% 
% _2.1 Permutation Testing on the regions_

% visualize bipolar relationship
bipolar_sig = {};
for i = 1:length(bipolar_timeseries)    
    bipolar_sig{i} = significant_correlations(bipolar_timeseries{i}, .05); %doing row here I changed the code
    figure;
    %swarmchart(ones(nnz(bipolar_sig{i}), 1), bipolar_sig{i}(bipolar_sig{i}~=0))
    %title("Bipolar Subject: " +i)
    %subtitle("Positive Values: "+nnz(bipolar_sig{i}>0) + " Negative Values: "+nnz(bipolar_sig{i}<0))
end

control_sig = {};
for i = 1:length(control_timeseries)    
    control_sig{i} = significant_correlations(control_timeseries{i}, .05); %doing row here I changed the code
    figure;
    %swarmchart(ones(nnz(control_sig{i}), 1), control_sig{i}(control_sig{i}~=0))
    %title("Control Subject: " +i)
    %subtitle("Positive Values: "+nnz(control_sig{i}>0) + " Negative Values: "+nnz(control_sig{i}<0))
end
%% 
% _2.2 Thresholding_
% 
% I have task data so I decided to permute the rows. The reason for this is 
% I want to have correlations that are signifiicant based on the time course of 
% the regions in which the task was administered. To say this another way I am 
% interested in the way that the regions of interest correlate with each other 
% over the time of the task. When I look at the data graphs that were generated 
% roughly a super majority of the significant values are positive with enough 
% data points i.e. signal to justify thresholding out the negative weights especially 
% if eventually I average out.

%threshold
for i = 1:length(bipolar_sig)
    bipolar_sig{i}(bipolar_sig{i}<0) = 0;
end

for i = 1:length(control_sig)
    control_sig{i}(control_sig{i}<0) = 0;
end
%% 
% *3. Network Analysis*
% 
% _3.1 Module Detection 
% 
% As a quick EDA step I would just like to visualize two of the random networks 
% before and after louvain clustering

randomnetwork = randi([1 14], 1);

vis_network(bipolar_sig{randomnetwork});
vis_network(control_sig{randomnetwork});
%% 
% The first step will be to run the community louvain algorithim to find an 
% intial paritipation amongst the control and patient groups.

bipolar_part = {}; 
gamma = 0.5:0.1:1.5;
for i = 1:14
    for g = 1:length(gamma)
        bipolar_part{i}(:,g) = consensus_community_louvain_with_finetuning(bipolar_sig{i}, gamma(g));
    end
end

control_part = {}; 
gamma = 0.5:0.1:1.5;
for i = 1:14
    for g = 1:length(gamma)
        control_part{i}(:,g) = consensus_community_louvain_with_finetuning(control_sig{i}, gamma(g));
    end
end

randompart = randi([1 length(gamma)], 1);
vis_network(bipolar_part{randomnetwork}(:,randompart));
vis_network(control_part{randomnetwork}(:,randompart));

tabulate(bipolar_part{randomnetwork}(:,randompart));
tabulate(control_part{randomnetwork}(:,randompart));
%% 
% _3.2 Network Metrics_ 
% 
% I am interested in the difference of networks of bipolar vs control patients. 
% Ultimately bipolar disorder is a disease that is very hard to diagnosis due 
% to a lack of quantitative tests. I would be interested in seeing if I can see 
% a difference in network metrics between groups. When examing a network using 
% various metrics I would like to look at this at two levels. The first will be 
% nodal level measures. They measures tell me something about how the nodes interact 
% on an individual basis. In this case the nodes will be the regions of interest 
% from the Schaffer 100 parcellation  atlas. So Im going to be interested in how 
% these regions interact with eachother across the board. I am also curious as 
% to how the network behaves overall. 
%% 
% First I will look at node level metrics specifically degree and participation 
% coefficient. I have choosen these metrics because they are simple yet powerful. 
% Degree is a measure that will allow me to examine how different regions interact 
% during a task while participation coefficent will tell me how important various 
% regions are to the overall integration of the network. When examing degree what 
% I am expecting to see will be that bipolar disorder will have similar degree 
% but lower participation coefficent. The reason for this is that I would guess 
% that bipolar nodes are connected to each other in a similar amount but not as 
% well as interconnected leading to the disease symptoms that we see. In other 
% words parts of the brain dont talk as much to each other as they should in a 
% healthly brain at a nodal level.  

%Strength
bipolar_str = [];
for i = 1:14
    bipolar_str = [bipolar_str, sum(bipolar_sig{i}, 2)];
end

control_str = [];
for i = 1:14
    control_str = [control_str, sum(control_sig{i}, 2)];
end

figure;
boxchart(bipolar_str);
title("Strength of Bipolar Nodes")
xlabel("Patients")

figure;
boxchart(control_str);
title("Strength of Control Nodes")
xlabel("Controls")

%Participation Coeff
bipolar_pc = {i};
for i = 1:14
    for q = 1:size(bipolar_part{i}, 2) 
        mod_deg = zeros(size(bipolar_part{i}, 1), max(bipolar_part{i}(:, q)));
        for p = 1:max(bipolar_part{i}(:, q))
            mod_deg(:, p) = sum(bipolar_sig{i}(:, bipolar_part{i}(:, q)==p), 2);
        end
        bipolar_pc{i} = 1 - sum((mod_deg ./ sum(mod_deg, 2)).^2, 2);
    end
end

control_pc = {};
for i = 1:14
    for q = 1:size(control_part{i}, 2) 
        mod_deg = zeros(size(control_part{i}, 1), max(control_part{i}(:, q)));
        for p = 1:max(control_part{i}(:, q))
            mod_deg(:, p) = sum(control_sig{i}(:, control_part{i}(:, q)==p), 2);
        end
        control_pc{i} = 1 - sum((mod_deg ./ sum(mod_deg, 2)).^2, 2);
    end
end

bipolar_pc = cell2mat(bipolar_pc);
control_pc = cell2mat(control_pc);

figure;
boxchart(bipolar_pc);
title("Participation Coefficient of Bipolar Nodes")
xlabel("Patients")

figure;
boxchart(control_pc);
title("Participation Coefficient of Control Nodes")
xlabel("Controls")
%% 
% Next I will look at network measures such as small worldness and rich club. 
% Ultimately these measures are much harder to interpert. However, small worldness 
% seems interesting in the context of task based data. I would be curious to see 
% how the network is structured when the brain is actually operating. Also, this 
% is a very important measure in the literature. I am expecting to see that bipolar 
% patients exhibit less small worldness. As I mentioned above I am expecting less 
% integration of the network across modules overall so it would make that the 
% bipolar patients a lower small world cofficient compared to the healthy patients. 
% Rich clubs given that they measure how other well connected nodes are connected 
% seem like a good complinteray measure to participation coefficent and degree. 
% I would get similar results between the two. In other words, I would expect 
% more signficant small worldness and higher rich club coefficents for the control 
% groups than the bipolar groups because I am expecting better integration of 
% the networks.

%Small Worldness
bipolar_transitivity = [];
bipolar_efficiency = [];
for i = 1:14    
    bipolar_transitivity = [bipolar_transitivity, transitivity_wu(bipolar_sig{i})];
    bipolar_efficiency = [bipolar_efficiency, efficiency_wei(bipolar_sig{i})];
end
bipolar_sw = bipolar_transitivity ./ bipolar_efficiency;

control_transitivity = [];
control_efficiency = [];
for i = 1:14    
    control_transitivity = [control_transitivity, transitivity_wu(control_sig{i})];
    control_efficiency = [control_efficiency, efficiency_wei(control_sig{i})];
end

control_sw = control_transitivity ./ control_efficiency;

figure;
scatter([1:14], bipolar_sw)
title("Bipolar Small Worldness")
subtitle("Count of small world networks: " + sum(bipolar_sw > 1))

xlabel('Patients')
ylabel('Small World Coefficient')

figure;
scatter([1:14], control_sw)
title("Control Small Worldness")
subtitle("Count of small world networks: " + sum(control_sw > 1))

xlim([0.0 14.0])
ylim([0.850 1.150])
xlabel('Control Groups')
ylabel('Small World Coefficient')

%Rich Club
bipolar_rich_club = [];
control_rich_club = [];

for i = 1:14
    r = rich_club_wu(bipolar_sig{i});
    bipolar_rich_club = [bipolar_rich_club, r(end)];
end

for i = 1:14
    r = rich_club_wu(control_sig{i});
    control_rich_club = [control_rich_club, r(end)];
end

bipolar_rich_club(isnan(bipolar_rich_club))=0;
control_rich_club(isnan(control_rich_club))=0;

figure;
scatter([1:14], bipolar_rich_club)
title("Bipolar Rich Clubs")
subtitle("Average Rich Club Coefficient: " + mean(bipolar_rich_club, "omitnan"))
xlim([1 14])
ylim([.5 1])
xlabel('Patients')
ylabel('Rich Club Coefficient')


figure;
scatter([1:14], control_rich_club)
title("Control Rich Clubs")
subtitle("Average Rich Club Coefficient: " + mean(control_rich_club, "omitnan"))
xlim([1 14])
ylim([.5 1])

xlim([1.00 14.00])
ylim([0.114 1.000])
xlabel('Control Group')
ylabel('Rich Club Coefficient')
%% 
% *4. Model Generation*
% 
% 1. Use simman to find models
% 
% 2. Generate a null distribution from each model
% 
% 2a. Generate a simman model of each person i.e. 28
% 
% 2b. For each model calculate a statistic of each person
% 
% 2c. Calculate average across group 
% 
% 2d. compare group averages
% 
% 2e. Repeat steps 2a - 2d N times 
% 
% 3. Compare empirical test statistic from group average differences 

samples = 10;

bipolar_model_1 = {};
bipolar_model_2 = {};

control_model_1 = {};
control_model_2 = {};

for i = 1:samples %every sample generates 14 of each
    bipolar_model_1{i} = simann(bipolar_timeseries, 1); 
    bipolar_model_2{i} = simann(bipolar_timeseries, 2);
    %bipolar_model_3(:,:,i) = simann(bipolar_timeseries, 3);
    %bipolar_model_4(:,:,i) = simann(bipolar_timeseries, 4);

    control_model_1{i} = simann(control_timeseries, 1);
    control_model_2{i} = simann(control_timeseries, 2);
    %control_model_3(:,:,i) = simann(control_timeseries, 3);
    %control_model_4(:,:,i) = simann(control_timeseries, 4);
end

bipolar_models = {bipolar_model_1, bipolar_model_2};
control_models = {control_model_1, control_model_2};
%%
null_distribution = {};
for i = 1:2
    bipolar_model = bipolar_models{i};
    control_model = control_models{i};
    current_null = [];
    for sample = 1:samples
        bipolar_sample = bipolar_model{sample};
        control_sample = control_model{sample};
    
        bipolar_richclub = 0;
        control_richclub = 0;
        
        for subject = 1:14

            bipolar_subject = bipolar_sample(:,:,subject);
            bipolar_thresh = threshold_proportional(bipolar_subject, .5);
            current_richclub = rich_club_wu(bipolar_thresh);
            current_richclub(isnan(current_richclub))=0;
            bipolar_richclub =+ current_richclub(end); 
    
            control_subject = control_sample(:,:,subject);
            control_thresh = threshold_proportional(control_subject, .5);
            current_richclub = rich_club_wu(control_thresh);
            current_richclub(isnan(current_richclub))=0;
            control_richclub =+ current_richclub(end);
        end

        current_null = [current_null; abs(bipolar_richclub - control_richclub) / 14];
    
    end

    null_distribution{i} = current_null;
end

for j = 1:14      
    bipolar_corr = corr(bipolar_timeseries{j}');
    bipolar_thresh = threshold_proportional(bipolar_corr, 0.5);
    current_richclub = rich_club_wu(bipolar_thresh);
    current_richclub(isnan(current_richclub))=0;
    bipolar_richclub =+ richclub(end); 
     
    control_corr = corr(control_timeseries{j}');
    control_thresh = threshold_proportional(control_corr, 0.5);
    current_richclub = rich_club_wu(control_thresh);
    current_richclub(isnan(current_richclub))=0;
    control_richclub =+ current_richclub(end);
end  
empirical_richclub = bipolar_richclub / 14 - control_richclub / 14;
%%
null_distribution = {};
for i = 1:2
    bipolar_model = bipolar_models{i};
    control_model = control_models{i};
    current_null = [];
    for sample = 1:samples
        bipolar_sample = bipolar_model{sample};
        control_sample = control_model{sample};
    
        bipolar_str = zeros(100, 1);
        control_str = zeros(100, 1);
        
        for subject = 1:14

            bipolar_subject = bipolar_sample(:,:,subject);
            bipolar_thresh = threshold_proportional(bipolar_subject, .5);
            strength = sum(bipolar_thresh, 2);
            bipolar_str =+ strength; 
    
            control_subject = control_sample(:,:,subject);
            control_thresh = threshold_proportional(control_subject, .5);
            strength = sum(control_thresh, 2);
            control_str =+ strength;
        end

        current_null = [current_null, abs(bipolar_str - control_str) / 14];
    
    end

    null_distribution{i} = current_null;
end

for j = 1:14
    bipolar_corr = corr(bipolar_timeseries{j}');
    bipolar_thresh = threshold_proportional(bipolar_corr, 0.5);
    strength = sum(bipolar_thresh, 2);
    bipolar_str =+ strength; 

    control_corr = corr(control_timeseries{j}');
    control_thresh = threshold_proportional(control_corr, .5);
    strength = sum(control_thresh, 2);
    control_str =+ strength;
end

empirical_strength = abs(bipolar_str  - control_str)  / 14;

%%
null_distribution = {};
for i = 1:2
    bipolar_model = bipolar_models{i};
    control_model = control_models{i};

    current_null = [];    
    for sample = 1:samples
        bipolar_sample = bipolar_model{sample};
        control_sample = control_model{sample};
        
        bipolar_pc_total = zeros(100, 1);
        control_pc_total = zeros(100, 1);
        for subject = 1:14
            gamma = 0.5:0.25:1; 
            bipolar_thresh= max(bipolar_sample(:,:,subject),0);
            control_thresh = max(control_sample(:,:,subject),0);
            for g = 1:length(gamma)
                bipolar_parts(:,g) = consensus_community_louvain_with_finetuning(bipolar_thresh, gamma(g));
                control_parts(:,g) = consensus_community_louvain_with_finetuning(control_thresh, gamma(g));              
            end
            part = randi(length(gamma));

            bipolar_part = bipolar_parts(:, part);        
            control_part = control_parts(:, part);

            for q = 1:size(bipolar_part, 2) 
                mod_deg = zeros(size(bipolar_part, 1), max(bipolar_part(:, q)));                    
                for p = 1:max(bipolar_part(:, q))
                    mod_deg(:, p) = sum(bipolar_thresh(:, bipolar_part(:, q)==p), 2);                
                end
                bipolar_pc = 1 - sum((mod_deg ./ sum(mod_deg, 2)).^2, 2);
            end

            for q = 1:size(control_part, 2) 
                mod_deg = zeros(size(control_part, 1), max(control_part(:, q)));                    
                for p = 1:max(control_part(:, q))
                    mod_deg(:, p) = sum(control_thresh(:, control_part(:, q)==p), 2);                
                end
                control_pc = 1 - sum((mod_deg ./ sum(mod_deg, 2)).^2, 2);
            end

            bipolar_pc_total =+ bipolar_pc;
            control_pc_total =+ control_pc;
        end 
        
        current_null = [current_null, abs(bipolar_pc_total - control_pc_total)/14];
    end 
    null_distribution{i} = current_null;
end
%% 
% *HELPER FUNCTIONS*
% 
% 

function [corr_data] = simann(data, constraint)  
    corr_data = zeros(100, 100, 14);
    for i = 1:length(data)
        if constraint == 1
            simann_ts = simann_model(data{i});
        elseif constraint == 2
            simann_ts = simann_model(data{i},'varnode', true, 'vartime', true);
        elseif constraint == 3
            simann_ts = simann_model(data{i}, 'varnode', true, 'vartime', true, 'covnodemode', true);
        elseif constraint == 4
            simann_ts = simann_model(data{i}, 'varnode', true, 'vartime', true, 'covnodemode', true,'covsystem', true);
        end 
        corr_data(:,:,i) = corr(simann_ts');
    end
end


function vis_timeseries(T)

% function to visualize timeseries and correlation matrices

figure,
tiledlayout(2, 1)

nexttile
imagesc(T)
xlabel('timepoints')
ylabel('brain regions')
colorbar

nexttile
imagesc(corr(T'), [-0.5 0.5])
xlabel('brain regions')
ylabel('brain regions')
axis square
colorbar
end

function vis_pvalues(C, P0, P1)
% function to visualize p values

figure,
tiledlayout(3, 1)

% visualize relationship with correlation
nexttile
plot(C(:), P0(:), '.'), hold on
if exist('P1', 'var')
    plot(C(:), P1(:), '.')
end
legend( ...
    'MATLAB p values', 'permutation p values', ...
    'location', 'bestoutside');
legend('boxoff')
xlabel('correlation')
ylabel('p value')
axis square
grid on

% visualize histograms
nexttile
histogram(P0, 0:0.05:1), hold on
if exist('P1', 'var')
    histogram(P1, 0:0.05:1),
end
legend( ...
    'MATLAB p values', 'permutation p values', ...
    'location', 'bestoutside');
legend('boxoff')
xlabel('p value')
axis square
grid on

% visualize relationship with each other
if exist('P1', 'var')
    nexttile,
    P0 = max(1e-10, P0);
    P1 = max(1e-10, P1);
    loglog(P0, P1, '.'); hold on
    loglog([1e-10 1], [1e-10 1], 'k')
    xlabel('MATLAB p values')
    ylabel('permutation p values')
    axis([1e-10 1 1e-10 1])
    axis square;
    grid on
end

end

function W = significant_correlations(T, pval)

% compute original correlation
[n, t] = size(T);
W = corr(T');

s = 1000;
P = zeros(n, n, s);
for i = 1:s
    T0 = T;
    
    % shuffle each row separately
    for j = 1:n
        T0(j, randperm(t)) = T(j, 1:t);
    end
    
    % compute null correlation;
    P(:, :, i) = (abs(corr(T0')) >= abs(W));
end
P = mean(P, 3);

% threshold by p value
W(P >= pval) = 0;

end

function [m1, q1_next] = community_louvain_with_finetuning(W, gamma)

    if ~exist('gamma', 'var')
        gamma = 1;
    end
    
    % number of nodes
    n = size(W, 1);
    
    % initial modularity vector
    m1 = 1:n;
    
    % initial modularity values
    q1_prev = -inf;
    q1_next = 0;
    
    %%% begin iterative finetuning %%%
    % while modularity increases
    while q1_next - q1_prev > 1e-5
        
        % run modularity with previous affiliation
        q1_prev = q1_next;
        [m1, q1_next] = community_louvain(W, gamma, m1);
    end
    %%% end iterative finetuning %%%
end


function [M, Q] = consensus_community_louvain_with_finetuning(W, gamma)

    if ~exist('gamma', 'var') % if gamma does not exist assign gamma to be one 
        gamma = 1;
    end
    
    n = size(W, 1); % mostly just shop keeping 
    k = 100;
    
    P = W;
    while 1
        M = zeros(n, k); % parition matrix
        Q = zeros(1, k);% modularity score matrix 
        for i = 1:k % over each column in M
            [M(:, i), Q(i)] = community_louvain_with_finetuning(P, gamma); % call the louvain algorithim to generate a consensus 
        end
        
        P = zeros(n); %set P to zero for a new set up
        for i = 1:k
            P = P + (M(:, i)==M(:, i)'); % calculate the consesus for P under 
        end
        P = P / k; %divide by k to have the max value be one i.e. a consensus 
        
        if all((P==0) | (P==1)) % if we have a consensus return the parition if not P will get passed back continuing the algorithim untill a consesus is reach
            break
        end
    end
    M = M(:, 1);
    Q = Q(1);

end

function vis_network(W)
    % function to network matrices and diagrams
    
    figure,
    tiledlayout(1, 2, 'padding', 'none')
    
    % visualize network matrix
    nexttile;
    imagesc(W);
    title('network matrix')
    axis square;
    colorbar;
    
    % visualize network diagram
    nexttile;
    A = max(W, 0);
    for i = 1:height(A)
        [~, ix] = sort(A(i, :), 'descend');
        A(i, ix(6:end)) = 0;
    end
    A = A + A';
    
    G = graph(A);
    plot(G, 'LineWidth', G.Edges.Weight/2)
    title('network diagram')
    axis square;

end