clc;
clear;
clf;
format long e;

small_A= load("small_ex.mat");
med_A1= load("medium_ex1.mat");
med_A2= load("medium_ex2.mat");

nnzmax = 100000000000;
tlimit = 1000000000000;

[L, flag, k] = spcholesky(small_A, nnzmax, tlimit);