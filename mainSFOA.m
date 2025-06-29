% main code 
clear all;close all;clc

% Parameters of algorithm
Npop = 50;              % Number of search agents
Max_it = 2000;          % Maximum number of iterations

func_num = 1;   % Function_number
Function_name = ['F',num2str(func_num)]; % Name of the test function
[lb,ub,nD,fobj]=Get_Functions_details(Function_name);

[xpos,fval,Curve] = SFOA(Npop,Max_it,lb,ub,nD,fobj); 