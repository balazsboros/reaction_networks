function out = ode_prlgrm_wide
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = @der4;
out{9} = @der5;

% --------------------------------------------------------------------------
function dydt = fun_eval(tau,kmrgd,par_a,par_b,par_t)
xstar=par_t;
ystar=par_t^2/4;
zstar=par_t^4/4;
z=(xstar+2*ystar+4*zstar-kmrgd(1)-2*kmrgd(2))/4;
q=1e2;
dydt=[q*2*kmrgd(2)*(z-kmrgd(1)^2*kmrgd(2));
q*(-kmrgd(2)*(z-kmrgd(1)^2*kmrgd(2))+2*(par_b*kmrgd(2)+par_a*kmrgd(1)^2)*(z-4*kmrgd(2)^2));];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(ode_prlgrm_wide);
y0=[0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(tau,kmrgd,par_a,par_b,par_t)
jac=[ -200*kmrgd(2)*(2*kmrgd(1)*kmrgd(2) + 1/4) , 50*par_t - 100*kmrgd(2) - 50*kmrgd(1) - 200*kmrgd(1)^2*kmrgd(2) + 25*par_t^2 + 50*par_t^4 - 200*kmrgd(2)*(kmrgd(1)^2 + 1/2) ; 100*kmrgd(2)*(2*kmrgd(1)*kmrgd(2) + 1/4) - 50*kmrgd(2)*par_b - 50*kmrgd(1)^2*par_a - 400*kmrgd(1)*par_a*(kmrgd(1)/4 + kmrgd(2)/2 - par_t/4 + 4*kmrgd(2)^2 - par_t^2/8 - par_t^4/4) , 25*kmrgd(1) + 50*kmrgd(2) - 25*par_t - 200*par_b*(kmrgd(1)/4 + kmrgd(2)/2 - par_t/4 + 4*kmrgd(2)^2 - par_t^2/8 - par_t^4/4) + 100*kmrgd(1)^2*kmrgd(2) - 100*(8*kmrgd(2) + 1/2)*(2*kmrgd(2)*par_b + 2*kmrgd(1)^2*par_a) - (25*par_t^2)/2 - 25*par_t^4 + 100*kmrgd(2)*(kmrgd(1)^2 + 1/2) ];
% --------------------------------------------------------------------------
function jacp = jacobianp(tau,kmrgd,par_a,par_b,par_t)
jacp=[ 0 , 0 , 200*kmrgd(2)*(par_t/4 + par_t^3 + 1/4) ; -200*kmrgd(1)^2*(kmrgd(1)/4 + kmrgd(2)/2 - par_t/4 + 4*kmrgd(2)^2 - par_t^2/8 - par_t^4/4) , -200*kmrgd(2)*(kmrgd(1)/4 + kmrgd(2)/2 - par_t/4 + 4*kmrgd(2)^2 - par_t^2/8 - par_t^4/4) , 100*(2*kmrgd(2)*par_b + 2*kmrgd(1)^2*par_a)*(par_t/4 + par_t^3 + 1/4) - 100*kmrgd(2)*(par_t/4 + par_t^3 + 1/4) ];
% --------------------------------------------------------------------------
function hess = hessians(tau,kmrgd,par_a,par_b,par_t)
hess1=[ -400*kmrgd(2)^2 , - 800*kmrgd(1)*kmrgd(2) - 50 ; 200*kmrgd(2)^2 - 200*kmrgd(1)*par_a - 400*par_a*(kmrgd(1)/4 + kmrgd(2)/2 - par_t/4 + 4*kmrgd(2)^2 - par_t^2/8 - par_t^4/4) , 400*kmrgd(1)*kmrgd(2) - 50*par_b - 400*kmrgd(1)*par_a*(8*kmrgd(2) + 1/2) + 25 ];
hess2=[ - 800*kmrgd(1)*kmrgd(2) - 50 , - 400*kmrgd(1)^2 - 200 ; 400*kmrgd(1)*kmrgd(2) - 50*par_b - 400*kmrgd(1)*par_a*(8*kmrgd(2) + 1/2) + 25 , 200*kmrgd(1)^2 - 1600*kmrgd(1)^2*par_a - 1600*kmrgd(2)*par_b - 400*par_b*(8*kmrgd(2) + 1/2) + 100 ];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
% --------------------------------------------------------------------------
function hessp = hessiansp(tau,kmrgd,par_a,par_b,par_t)
hessp1=[ 0 , 0 ; - 50*kmrgd(1)^2 - 400*kmrgd(1)*(kmrgd(1)/4 + kmrgd(2)/2 - par_t/4 + 4*kmrgd(2)^2 - par_t^2/8 - par_t^4/4) , -200*kmrgd(1)^2*(8*kmrgd(2) + 1/2) ];
hessp2=[ 0 , 0 ; -50*kmrgd(2) , 50*par_t - 100*kmrgd(2) - 50*kmrgd(1) - 800*kmrgd(2)^2 + 25*par_t^2 + 50*par_t^4 - 200*kmrgd(2)*(8*kmrgd(2) + 1/2) ];
hessp3=[ 0 , 50*par_t + 200*par_t^3 + 50 ; 400*kmrgd(1)*par_a*(par_t/4 + par_t^3 + 1/4) , 200*par_b*(par_t/4 + par_t^3 + 1/4) - 25*par_t - 100*par_t^3 - 25 ];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
hessp(:,:,3) =hessp3;
%---------------------------------------------------------------------------
function tens3  = der3(tau,kmrgd,par_a,par_b,par_t)
tens31=[ 0 , -800*kmrgd(2) ; -300*par_a , 400*kmrgd(2) - 400*par_a*(8*kmrgd(2) + 1/2) ];
tens32=[ -800*kmrgd(2) , -800*kmrgd(1) ; 400*kmrgd(2) - 400*par_a*(8*kmrgd(2) + 1/2) , 400*kmrgd(1) - 3200*kmrgd(1)*par_a ];
tens33=[ -800*kmrgd(2) , -800*kmrgd(1) ; 400*kmrgd(2) - 400*par_a*(8*kmrgd(2) + 1/2) , 400*kmrgd(1) - 3200*kmrgd(1)*par_a ];
tens34=[ -800*kmrgd(1) , 0 ; 400*kmrgd(1) - 3200*kmrgd(1)*par_a , -4800*par_b ];
tens3(:,:,1,1) =tens31;
tens3(:,:,1,2) =tens32;
tens3(:,:,2,1) =tens33;
tens3(:,:,2,2) =tens34;
%---------------------------------------------------------------------------
function tens4  = der4(tau,kmrgd,par_a,par_b,par_t)
tens41=[ 0 , 0 ; 0 , 0 ];
tens42=[ 0 , -800 ; 0 , 400 - 3200*par_a ];
tens43=[ 0 , -800 ; 0 , 400 - 3200*par_a ];
tens44=[ -800 , 0 ; 400 - 3200*par_a , 0 ];
tens45=[ 0 , -800 ; 0 , 400 - 3200*par_a ];
tens46=[ -800 , 0 ; 400 - 3200*par_a , 0 ];
tens47=[ -800 , 0 ; 400 - 3200*par_a , 0 ];
tens48=[ 0 , 0 ; 0 , 0 ];
tens4(:,:,1,1,1) =tens41;
tens4(:,:,1,1,2) =tens42;
tens4(:,:,1,2,1) =tens43;
tens4(:,:,1,2,2) =tens44;
tens4(:,:,2,1,1) =tens45;
tens4(:,:,2,1,2) =tens46;
tens4(:,:,2,2,1) =tens47;
tens4(:,:,2,2,2) =tens48;
%---------------------------------------------------------------------------
function tens5  = der5(tau,kmrgd,par_a,par_b,par_t)
tens51=[ 0 , 0 ; 0 , 0 ];
tens52=[ 0 , 0 ; 0 , 0 ];
tens53=[ 0 , 0 ; 0 , 0 ];
tens54=[ 0 , 0 ; 0 , 0 ];
tens55=[ 0 , 0 ; 0 , 0 ];
tens56=[ 0 , 0 ; 0 , 0 ];
tens57=[ 0 , 0 ; 0 , 0 ];
tens58=[ 0 , 0 ; 0 , 0 ];
tens59=[ 0 , 0 ; 0 , 0 ];
tens510=[ 0 , 0 ; 0 , 0 ];
tens511=[ 0 , 0 ; 0 , 0 ];
tens512=[ 0 , 0 ; 0 , 0 ];
tens513=[ 0 , 0 ; 0 , 0 ];
tens514=[ 0 , 0 ; 0 , 0 ];
tens515=[ 0 , 0 ; 0 , 0 ];
tens516=[ 0 , 0 ; 0 , 0 ];
tens5(:,:,1,1,1,1) =tens51;
tens5(:,:,1,1,1,2) =tens52;
tens5(:,:,1,1,2,1) =tens53;
tens5(:,:,1,1,2,2) =tens54;
tens5(:,:,1,2,1,1) =tens55;
tens5(:,:,1,2,1,2) =tens56;
tens5(:,:,1,2,2,1) =tens57;
tens5(:,:,1,2,2,2) =tens58;
tens5(:,:,2,1,1,1) =tens59;
tens5(:,:,2,1,1,2) =tens510;
tens5(:,:,2,1,2,1) =tens511;
tens5(:,:,2,1,2,2) =tens512;
tens5(:,:,2,2,1,1) =tens513;
tens5(:,:,2,2,1,2) =tens514;
tens5(:,:,2,2,2,1) =tens515;
tens5(:,:,2,2,2,2) =tens516;
