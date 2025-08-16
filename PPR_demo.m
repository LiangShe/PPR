
clear

%% simulate model neuron

load stimuli.mat
stim = stim(9:24,9:24,1:5000); % use smaller stimuli for faster demo
FiringRate=50;
FrameRate=20;
[resp, rf1, rf2] = Model_Neuron_V1( stim, 'complex', FiringRate, FrameRate ); 
resp = resp';
fprintf('total spikes %f\n',sum(resp));

%% compute PPR

outputpath = 'results';

% weight for the penalize term of nonsmoothness 
lamda = [0  1e-7  5e-7  1e-6  5e-6  1e-5  5e-5 ];

for i=1:length(lamda)
    fprintf("lamda = %e\n", lamda(i))
    PPR_2order_4(stim, resp, lamda(i), 'complex_model', outputpath); close all
end

%% view results
% TODO
