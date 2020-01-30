clc;
clear;
clf;
format long e;

Tm = [2 -1  0  0  0;
      -1 2 -1  0  0;
      0 -1  2 -1  0;
      0 0  -1  2 -1;
      0 0   0 -1  2;]
  I = eye(5);
  D = Tm+2*I;
  F = -I;
  Z = zeros(5);
  
  Tmm = [D F Z Z Z;
         F D F Z Z;
         Z F D F Z;
         Z Z F D F;
         Z Z Z F D;]
     
  
  