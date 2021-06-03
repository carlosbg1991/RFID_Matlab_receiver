function y=tsallis(x,p,w,q)
% Tsallis distribution
% "Tsallis model-based separation of overlapped peak signals", Li YuanLu,
% Zhang YingChao, and Tang HuiQiang, Science Shina, 53, 823-832, 2010.
% http://link.springer.com/article/10.1007/s11432-010-0054-4#page-1
% where p=mu and w=sigma
% 1<q<3.   Shape->Gaussian as q->1; Shape->Lorentzian as q->2
y=(1+((q-1)/(3-q)).*(((x-p).^2)./(w.^2))).^(-1./(q-1));