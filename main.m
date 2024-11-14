clc;
clear;

datafiles = {'mul_dermatology'};

for datafile = datafiles

    load([datafile{1} '.mat']);
    output_file = ['results/' datafile{1} '_results.mat'];
    fprintf('Dataset:%s\t',datafile{1});

    addpath(genpath('./'));


    nv=length(data); % view
    ns=length(unique(labels)); % k
    for ni=1:nv
        data{ni,1} = mapstd(data{ni}',0,1); % mapminmax(data{ni}',0,1);%Normalize X
    end

    anchor_traverse = 4;
    numanchor = (1:anchor_traverse)*ns;
    alpha = 0.1;beta = 1;lan = 1;

    options.alpha = alpha;
    options.beta = beta;
    options.lan = lan;

    [finU,alfa,time,obj] = AMVFC_FDAG(data,options,numanchor, ns);
    [pred_labels, ~] = litekmeans(finU,ns,'MaxIter', 100,'Replicates',10);

    result_cluster = ClusteringMeasure(labels, pred_labels);
    nmi = result_cluster(2);
    acc = result_cluster(1);
    purity = result_cluster(3);
    ARI = result_cluster(4);


end