% Funktion til at beregne iskoncentrationer
% Metode: Nasa Team Sea Ice Algorithm
% Input: TA(6V, 6H, 10V, 10H, 18V, 18H, 23V, 23H, 37V, 37H)
% Output: C i %
%***********************************
% Dorthe Hofman-Bang, 19. september 2003
%***********************************
function[C_FY, C_MY, CT]=is_C_funktion_AMSR(TA)

TB19H=TA(:,6);      %load brightness temperatures for algorithm
TB19V=TA(:,5);
TB37V=TA(:,9);

%************************************
%konstanter

a0=3286.56;             
b0=-790.321;
c0=2032.20;

a1=-20764.9;
b1=13825.8;
c1=9241.50;

a2=23893.1;
b2=-33104.7;
c2=-5655.62;

a3=47944.5;
b3=-47720.8;
c3=-12864.9;

%*******************************
%NASA Team Algorithm

TPR=TB19V-TB19H;
NPR=TB19V+TB19H;

PR=TPR./NPR;

TGR=TB37V-TB19V;
NGR=TB37V+TB19V;

GR=TGR./NGR;

D=c0+c1.*PR+c2.*GR+c3.*PR.*GR;

C_FY=((a0+a1.*PR+a2.*GR+a3.*PR.*GR)./D);  %First Year concentration
C_MY=((b0+b1.*PR+b2.*GR+b3.*PR.*GR)./D);  %Multi Year concentration

CT=(C_MY+C_FY); %The total iceconcentration 
