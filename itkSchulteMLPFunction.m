function [] = itkSchulteMLPFunction()
    clc;
    figure(1);
    CompBetaPSquare();
    figure(2);
    CompareSpatialSigma();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compares the fits for the function 1/beta²p². See [Schulte, 2008].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = CompBetaPSquare()

    u=[0:0.1:20];

    % Schulte
    a=zeros(6,1);
    a(1)=7.457e-6;
    a(2)=4.548e-7;
    a(3)=-5.777e-8;
    a(4)=1.301e-8;
    a(5)=-9.228e-10;
    a(6)=2.687e-11;
    schulte=zeros(1,length(u));
    for i=0:5
        schulte = schulte + a(i+1) .* u.^i;
    end

    % Williams
    a(1)=7.507e-4;
    a(2)=3.320e-5;
    a(3)=-4.171e-7;
    a(4)=4.488e-7;
    a(5)=-3.739e-8;
    a(6)=1.455e-9;
    a=a.*0.01;
    williams=zeros(1,length(u));
    for i=0:5
        williams = williams + a(i+1) .* u.^i;
    end

    % Rit
    a(1)=7.444724e-06;
    a(2)=5.463937e-07;
    a(3)=-9.986645e-08;
    a(4)=2.026409e-08;
    a(5)=-1.420501e-09;
    a(6)=3.899100e-11;
    rit=zeros(1,length(u));
    for i=0:5
        rit = rit + a(i+1) .* u.^i;
    end

    clf;
    hold on;
    plot(u,schulte, 'r-');
    plot(u,williams, 'b-');
    plot(u,rit, 'k-');
    % 
    % E0=13.6;
    % X0=36.1;
    % sigmasq = zeros(1,length(u));
    % for i=0:5
    %     sigmasq = sigmasq + a(i+1)./(i+1) .* u.^(i+1);
    % end
    % sigmasq = sigmasq .* E0*E0.*(1+0.038*log(u/X0))./X0;
    % 
    % plot(u, sqrt(sigmasq), 'r-');
    
function [] = CompareSpatialSigma()
    % Computed with code
    calc = [0. 0.0316514 0.0942212 0.179374 0.284538 0.408824 0.552243 0.71532 0.898855 1.10379 1.33115 1.58214 1.85821 2.16138 2.49446 2.86144 3.26793 3.72148 4.23201 4.81203 5.47676];
    % Computed with MC simulation without hadronic processes
    sim_wo_had = [0.0620 0.0538 0.0969 0.1710 0.2793 0.3695 0.4901 0.6299 0.7878 0.9545 1.1308 1.3198 1.5202 1.7364 1.9602 2.1954 2.4443 2.7032 2.9780 3.2621 3.5602];
    % Computed with MC simulation using all protons
    sim = [0.4890 0.7338 1.2242 1.7920 2.3543 2.9185 3.4722 3.9982 4.5000 4.9619 5.3770 5.7720 6.1082 6.4049 6.6641 6.8830 7.0589 7.2042 7.3460 7.4762 7.6156];
    % Idem without secondaries
    sim_wo_sec = [0.0614 0.1581 0.3569 0.5958 0.8563 1.1270 1.4039 1.6906 1.9862 2.2897 2.5959 2.9143 3.2370 3.5669 3.9037 4.2510 4.6103 4.9829 5.3770 5.7924 6.2356];
    clf;
    hold on;
    u=[0:10:200];
    plot(u,calc, 'r-');
    plot(u,sim, 'b-');
    plot(u,sim_wo_had, 'g-');
    plot(u,sim_wo_sec, 'k-');



% plot(u,(calc./sim), 'r-');
polyfit(u,(calc./sim),1)