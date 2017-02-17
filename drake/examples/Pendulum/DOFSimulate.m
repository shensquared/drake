function DOFSimulate()

    pv = PendulumVisualizer();
    sys=MyPendulumCL;
    
    initial_state=[0;0]+0.2*randn(2,1);
    % if the estimation and the initial state are equal at nonzero values, then it's possible for the both of them to deviate from the 
    % equlibrim at first and then start to converge again to the zero values.
    initial_estimate=initial_state;
    co_variance = .2;
    % initial_estimate(2)=initial_state(2)+covariance*randn(1,1);
    initial_augmented=[initial_state;initial_estimate];
    [ytraj,xtraj]=simulate(sys,[0,12],initial_augmented);
%     pv.playback(ytraj); 
    disp('done')
    figure(1);
    subplot(2,2,1);
    fnplt(xtraj,1);
    % figure(2);
    subplot(2,2,2);
    fnplt(xtraj,2);
    % figure(3);
    subplot(2,2,3);
    fnplt(xtraj,3);
    % figure(4);
    subplot(2,2,4);
    fnplt(xtraj,4);

end