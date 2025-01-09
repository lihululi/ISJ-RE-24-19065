clc
clear
close all

tic
M = 40;
num_chrom = 1;
N_chrom = 26;
Miter = 100;
Pm = 0.2;
Pc = 0.35;
Er = 0.05;
visualization = 1;

[BestChrom, cgcurve, fitness_ave] = Eg1_GA(M, num_chrom, Miter, Pc, Pm, Er, visualization);

disp('The best chromosome found: ')
BestChrom.Gene
disp('The best fitness value: ')
BestChrom.Fitness

R = BestChrom.R;
gene = BestChrom.Gene;
toc
