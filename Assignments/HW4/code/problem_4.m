set(0,'DefaultFigureWindowStyle','docked');
clc;
clear;
close all;
clf;
format long e;


%a)
load('HW4_P4a.mat');

T = nonsymmetric_lanczos(A, r, c, 20);
norm(eig(T) - eig(A))
