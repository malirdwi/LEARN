addpath('..');
clear all;   

Gamma =[
-1 1 0 0 0 0 0 1;
-1 1 0 0 0 0 1 0;
 1 -1 -1 0 0 0 0 0;
 0 0 0 -1 1 0 1 0;
 0 0 0 -1 1 0 0 1;
 0 0 0 1 -1 -1 0 0;
 0 0 1 0 0 0 -1 0;
 0 0 0 0 0 1 0 -1];
LEARNmain(Gamma)
