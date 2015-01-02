# Adaptive Control of Robot Manipulator

This repository contains matlab/simulink codes for my article "Adaptive Control of 4-DoF Robot manipulator". More precise, there you will a model of 4-DoF robot manipulator, codes for derriving dynamics of such manipulator symbollically and simulation model where robots is controllled by adaptive torque controller, which estimates some parametes (CoM of sub-bodies) of robots online.

Depends on : Peters Corke's Robotics Toolbox @ http://www.petercorke.com

projest.m  
Main file - modelling robot kinematics and dynamics, derrivation of dynamics symbollically

W_matrix.m
D_matrix.m
matrixes produced by project.m

parametrize.m - symbollical parameters linearization of W matrix 

W_adapt_matrix.m - torque controller with adaptive estimation of parameters

sl_robot.slx - simulink simulation of manipulator control with parameter estimation

write_fcn.m

Pavel Mironchyk
