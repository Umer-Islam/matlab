a = 5;
b = 10;

% division
a/b
%left division, it changes the values
a\b

% example 1-1a

 a = 5*(3/4) + (9/5)

 % example 1-1b

 (4^3)* ((3/4)+(9/(2*3)))

 %Assignment Operator "="
 x  = 34-6;

 %suppressing output-> use semiColon at the 
 % end to suppess the output on command window
x =2; y=4;z = x*y

%use $ who to see all the variables that we have


% by $whos will give us variables,
% their values,type and size(as vector/matrix)



%use the threeDots "..."(ellipsis) to extend a long line
FirstClassHolders = 72;
 Coach = 121;
 Crew = 8;
 TotalPeopleOnPlane = FirstClassHolders + Coach...
+ Crew

 %testing short vs long format
 format long
format_long = 3 + 11/26 + 2^1.2
 format short
  format_short = 3 + 11/26 + 2^1.2
 %we we dont want to round off the value
 format bank
 format_bank = 3 + 11/26 + 2^1.2


 % when we don't want to use exponent.

 no_short_e = 7.2 *3.2
 format short e
 format_short_e = 7.2 *3.2

 % To find the ratioanl expression or 😲 closest rational expression
 format rat
 value = .6
 

%Example 1-2: find the volume of sphere of radus 2 meter.
% where v = (4/3)*pi*r^3
radius_sphere = 2;
format long % to convert fraction value to decimal
volume_sphere = (4/3) * (pi) * radius_sphere^3



%sqrt of a value

square_root = sqrt(16)

%natural log of a value
format short
log(3.2)

% to find log with a base
log10(3)

%basic trig functions
cos(pi/4)

%use a before the trig function to use inverse like: acos, atan
format rat
atan(pi/3)

%complex numbers
complex_a = 2+3i;
complex_b = 1-i;
format default
complex_sum = complex_a + complex_b


 format rat
sin_rat = sin(pi/3)