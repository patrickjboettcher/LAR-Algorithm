function out=fdn_biquad(in, A, D, HB, HA)

% Feedback Delay Network
% 
% y=fdn_biquad(in, A, D, HB, HA) performs a feedback delay network on input
%   signals in. The number of delay paths, Dcnt, is determined by the length 
%   of D. Accordingly, Dcnt input signals are expected and Dcnt output signals
%   are returned. Additionally, fdn_biquad performs a biquad filtering on 
%   each of the delay path according to the filter coefficients HB and HA, 
%   both corresponding to B and A from the function filter. 
% 
% in: input signals, [length x Dcnt]
% A: feedback matrix, [Dcnt x Dcnt]
% D: delay vector (samples), [Dcnt x 1]
% HB: non-recursive coefficients of the biquad filters, [3 x Dcnt]
% HA: recursive coefficients of the biquad filters, [3 x Dcnt], see filter for details
% 
% out: output signals [length x Dcnt]
%
% 2021 Piotr Majdak, under EUPL v1.2

Dcnt=length(D);    % # of delay paths | was max(D) before
Dmax=max(D)+2;
N=length(in)+Dmax;
y=zeros(N,Dcnt);
s=zeros(N,Dcnt);
d=zeros(N,Dcnt);
x=zeros(N,Dcnt);
x(Dmax+1:N,:)=in; % was x(Dmax+1:N,:)=in;

for n=Dmax+1:N
  for dd=1:Dcnt
    s(n,dd)=sum(A(dd,:)*y(n-1,:)')+x(n,dd);
    d(n,dd)=s(n-D(dd),dd);
    y(n,dd)=HB(1,dd)*d(n,dd)+HB(2,dd)*d(n-1,dd)+HB(3,dd)*d(n-2,dd)-HA(2,dd)*y(n-1,dd)-HA(3,dd)*y(n-2,dd);
  end
end

out=y(Dmax+1:N,:);
