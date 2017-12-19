function [ delay ] = calc_wire_delay( tech, vdd_real, wire_level )

%switching occurs at 0.5 VDD
a=0.4;
b=0.7;

% tech=45;
% vdd_real=1.0;
verbose=0;
% wire_level=2; %1: M1, 2: intermediate 3: global

[ out_res, in_cap, par_cap, wire_res, wire_cap ] = calc_cacti_parameters( tech, vdd_real, verbose );

r=wire_res(wire_level)*1e6; %scale ohm/micron to ohm/meter
c=wire_cap(wire_level)*1e6; %scale F/micron to F/meter
r0=out_res;
c0=in_cap;
cp=par_cap;

%130nm from past literature (UT  Austin and Univ. Illinois Urbana)
% r=61*1e-3/1e-6;
% c=0.359*1e-15/1e-6;
% r0=8.8*1e3;
% c0=0.94*1e-15;
% cp=2.29*1e-15;

% cp=c0;

lopt=sqrt((b/a)*r0*c0*(1+(cp/c0)))/sqrt(r*c);
t=2*sqrt(r*c*r0*c0)*(b+sqrt(a*b*(1+(cp/c0))));

% display(['lopt = ' num2str(lopt*1e3) ' mm'])
%coskun paper says delay is 183 ps/mm
% display(['delay = ' num2str(t*1e12/1e3) ' ps/mm'])

delay=t*1e12/1e3; %ps/mm

end