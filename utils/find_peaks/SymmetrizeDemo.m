  % SymmetrizeDemo runs the 5 examples in the symmetrize.m help file
  disp('Example 1: Symmetrization of an exponentially broadend Gaussian.')
  t=[0:.01:20]';
  tau=1; % tau is exponential decay constant in x units (1 segment)
  y=expgaussian2(t,8,1,tau)'; 
  factor=1;% For Gaussian, factor is equal to tau
  symy=symmetrize(t,y,1); % No smoothing because signal is noiseless
  figure(1);clf;plot(t,y,t,symy)
  OriginalArea=trapz(t,y)
  SharpenedArea=trapz(t,symy)
  xlabel('t')
  ylabel('y')
  title('Blue: Original Gaussian    Red: Symmetrized by first derivative addition') 
 
    disp('Example 2: Symmetrization of an noisy exponentially broadended Lorentzian.')
  increment=5;t=0:increment:3000;
  tau=50; % tau is exponential decay constant in x units (1 segment)
  y=explorentzian(t,1000,100,-tau)'+.002.*randn(size(t)); % Position=1000, width=100
  factor=250; % For lorentzian, factor is about increment*tau
  symy=symmetrize(t,y,250,21,3); % smoothwidth=21, type=3
  figure(2);clf;plot(t,y,t,symy)
  xlabel('x')
  ylabel('y')
  title('Blue: Original signal    Red: Symmetrized by first derivative addition') 
  OriginalArea=trapz(t,y)
  SharpenedArea=trapz(t,symy)
  xlabel('t')
  ylabel('y')
  title('Blue: Original Lorentzian    Red: Symmetrized by first derivative addition') 
 
  disp('Example 3: Symmetrization of a signal with multiple EMG peaks')
  load DataMatrix3; % Data file from the SPECTRUM.ZIP archive
  x=DataMatrix3(:,1)';
  y=DataMatrix3(:,2)'; 
  tau=17; % tau is exponential decay constant in x units (1 segment)
  smoothwidth=5;
  symy=symmetrize(x,y,tau,smoothwidth); 
  figure(3);clf;plot(x,y,x,symy)
  xlabel('x')
  ylabel('y')
  title('Blue: Original signal    Red: Symmetrized by first derivative addition') 
 
  disp('Example 4: 8-segment symmetrization of 5 Gaussians with different exponential broadening time constants.')
  x=5:.1:65;
  y=modelpeaks2(x, [1 5 5 5 5], [1 1 1 1 1], [10 20 30 40 50], [3 3 3 3 3], [0 -5 -10 -15 -20]);
  factor=[.4 .5 .6 1.1 1.3 1.5 2 2]; % Defines 8-segment symmetrization factor
  Symmy=symmetrize(x,y,factor,3,1,1);
  figure(4);clf;plot(x,y,x,Symmy)
  xlabel('x')
  ylabel('y')
  title('Blue: Original signal    Red: Symmetrized by first derivative addition') 
 
  disp('Example 5: 30-segment gradient symmetrization of 5 Gaussians with different exponential broadening time constants.')
  x=5:.1:65;
  y=modelpeaks2(x, [1 5 5 5 5], [1 1 1 1 1], [10 20 30 40 50], [3 3 3 3 3], [0 -5 -10 -15 -20]);
  NumSegments=30; % Number of segments
  startfactor=0;
  endfactor=2.6;
  factorstep=(endfactor-startfactor)/NumSegments;
  factor=startfactor:factorstep:endfactor; % Defines multiple-segment symmetrization factor
  Symmy=symmetrize(x,y,factor,3,1,1);
  figure(5);clf;plot(x,y,x,Symmy)
  xlabel('x')
  ylabel('y')
  title('Blue: Original signal    Red: Symmetrized by first derivative addition')