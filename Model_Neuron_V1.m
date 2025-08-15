function [spikes, gab1, gab2, P] = Model_Neuron_V1( stim, neuron_type, FiringRate, FrameRate )
% Receptive field model neuron for V1 simple and complex cell
% 

% defalut parameters
if nargin<1
    stim=randn(16, 16, 10000);  % defalut stimulus White noise
end

if nargin<2
    neuron_type = 'simple';
end

if nargin<3
    FiringRate = 50;  % defalut FiringRate Hz
end
 
if nargin<4
    FrameRate = 20;  % defalut FrameRate Hz
end


% nonlinear
nf = 2; %  NonlinearFactor
threshold = 0.000; % threshold

[X, Y, Z]=size(stim);
time = Z/FrameRate;  % sec

% adaptaion to absolute intensity
% for those have similar average intensity throughout the stimuli
stim=stim-mean(stim(:)); 

% for those movies whose average intensity vary throughout(e.g. some natural stim)
% more precisely, average of stim(t-t0:t)should subtract from the stim(t)
 

%% temporal response function
a=4;  % ms
b=8;  % ms
t=5:1000/FrameRate:150;
trf=t/(a^2).*exp(-t/a)-(t/b^2).*exp(-t/b);  

%% spatial RF gabor
sigma=0.1;   % rf size
k=0.25;      % spatial frequency percentage/cycle
phi=0;       % spatial phase
theta=30;     % orientation, clockwise
asp=0.5;     % aspect ratio, determines ellipticity.
gab1=Gabor(sigma, k, phi, theta, asp, Y, X);
gab2=Gabor(sigma, k, phi+90, theta, asp, Y, X);


%% spatial temporal
strf1=zeros(X,Y,length(trf));
strf2=zeros(X,Y,length(trf));
for i=1:X
    for j=1:Y
        strf1(i,j,:)=gab1(i,j)*trf;
        strf2(i,j,:)=gab2(i,j)*trf;
    end
end


%% response
switch neuron_type
    
    case 'simple'
        resp1=zeros(1,Z+length(trf)-1);
        st=zeros(1,Z);
        f1=zeros(1,length(trf));
        for i=1:X
            for j=1:Y
                st(:)=stim(i,j,:);
                f1(:)=strf1(i,j,:);
                resp1=resp1+conv(st,f1);
            end
        end
        resp1 = max(resp1, threshold);
        resp = resp1(1:Z).^2;
        
    case 'complex'
        resp1=zeros(1,Z+length(trf)-1);
        resp2=zeros(1,Z+length(trf)-1);
        st=zeros(1,Z);
        f1=zeros(1,length(trf));
        f2=zeros(1,length(trf));
        for i=1:X
            for j=1:Y
                st(:)=stim(i,j,:);
                f1(:)=strf1(i,j,:);
                f2(:)=strf2(i,j,:);
                resp1=resp1+conv(st,f1);
                resp2=resp2+conv(st,f2);
            end
        end
        
        resp = resp1(1:Z).^2+resp2(1:Z).^2;

end
 
%% Poisson generator
P = resp/sum(resp)*FiringRate*time ; % normalize,

spikes = poissrnd(P);


 
%% Gabor fuction
function [ gabor ] = Gabor( sigma, k, phi, theta, asp, X, Y )
%INPUT
% sigma  % rf size
% k      % spatial frequency percentage/cycle
% phi    % spatial phase /degree
% theta  % orientation, clockwise
% asp    % aspect ratio, determines ellipticity.
% X  resolution
% Y
if nargin<1
sigma=0.1;   % rf size
end
if nargin<2
k=0.25;     % spatial frequency, percentage of the whole RF/cycle
end
if nargin<3
    phi=0;    % spatial phase /degree
end
if nargin<4
    theta=0;  % orientation, clockwise
end
if nargin<5
    asp=0.5;   % aspect ratio, determines ellipticity.
end
 
theta=theta*pi/180;
phi=phi*pi/180;
[xx,yy] = meshgrid(1:X,1:Y);
xx=xx/X; % normalize to 1
yy=yy/X; % to make the x/y ratio fixed 
yy=yy(end:-1:1,:);
 
x=(xx-0.5).*cos(theta)-(yy-0.5).*sin(theta);
y=(xx-0.5).*sin(theta)+(yy-0.5).*cos(theta);
gabor=1/2*pi/ sigma^2 /asp * exp(( -x.^2 -asp^2 * y.^2)/(2*sigma^2)).* cos(2*pi/k*x - phi);
 
