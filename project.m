%% Adaptive Control Project
% 
% Pavel Mironchyk 2014

slCharacterEncoding('ISO-8859-1');

clear all;

% parameters
L_upper_leg=0.5;
L_lower_leg=0.5;
%L_upper_leg
%L_lower_leg

%s='Rz(q_hip_z).Rx(q_hip_x).Ry(q_hip_z).Tz(L_upper_leg).Ry(q_knee_y).Tz(L_lower_leg).Ry(q_ankle_y).Tx(L_feet)';

% create the leg complex_leg_links based on DH parameters
%                               theta      d    a               alpha  
complex_leg_links(1)   = Link([    0       0    0.01             pi/2 ], 'standard');
complex_leg_links(1).m = 0.0;%0.1;
complex_leg_links(1).r = [0 0 0];
complex_leg_links(1).I = zeros(3,3);
complex_leg_links(1).G = 1000;
complex_leg_links(1).Jm = 200e-6;
complex_leg_links(1).B = 0;

complex_leg_links(2) = Link([    0       0    0.01             pi/2 ], 'standard');
complex_leg_links(2).m = 0.1;
complex_leg_links(2).r = [-0.005 0 0];
complex_leg_links(2).I = zeros(3);%eye(3)*2.0/5.0*0.1*0.0001;
complex_leg_links(2).G = 1000;
complex_leg_links(2).Jm = 200e-6;
complex_leg_links(2).B = 0;


complex_leg_links(3) = Link([    0       0    L_upper_leg     pi/2 ], 'standard');
complex_leg_links(3).m = 1;
complex_leg_links(3).r = [-L_upper_leg/2 0 0];
complex_leg_links(3).I = zeros(3);%eye(3)*2.0/5.0*1*0.001;
complex_leg_links(3).G = 1000;
complex_leg_links(3).Jm = 200e-6;
complex_leg_links(3).B = 0;

complex_leg_links(4) = Link([    0       0    L_lower_leg     0   ], 'standard');
complex_leg_links(4).m = 1;
complex_leg_links(4).r = [-L_lower_leg/2 0 0];
complex_leg_links(4).I = zeros(3); %eye(3)*2.0/5.0*1*0.001;
complex_leg_links(4).G = 1000;
complex_leg_links(4).Jm = 200e-6;
complex_leg_links(4).B = 0;



complex_leg_offsets = [0,  pi/2, 0, 0];

% now create a robot to represent a single leg
complex_leg = SerialLink(complex_leg_links, 'name', 'complex_leg', 'offset', complex_leg_offsets);

% q_hip_z, q_hip_y, q_hip_x, q_knee_y
%complex_leg.plot([0,0,-pi/10,pi/5]);
%set(gca, 'Zdir','reverse'); 

qz=[0,0,0,0];
complex_leg.accel([1.0, 10, 0, 0],[10.0, 10, 0, 0],qz)

% computing dynamics symbolically
syms q_hip_z q_hip_y q_hip_x q_knee_y real
syms dq_hip_z dq_hip_y dq_hip_x dq_knee_y real
syms ddq_hip_z ddq_hip_y ddq_hip_x ddq_knee_y real
syms l_upper_leg l_lower_leg real
syms mass_upper_leg mass_lower_leg real
syms com_upper_leg com_lower_leg real
syms g real

q  = [q_hip_z,  q_hip_y,  q_hip_x,  q_knee_y]';
dq = [dq_hip_z, dq_hip_y, dq_hip_x, dq_knee_y]';
ddq = [ddq_hip_z, ddq_hip_y, ddq_hip_x, ddq_knee_y]';


%%derriviation of Euler-Langrange dynamics
r_upper_body=rotz(q_hip_z)*roty(-q_hip_y)*rotx(q_hip_x)*[0; 0; com_upper_leg];
r_lower_body=rotz(q_hip_z)*roty(-q_hip_y)*rotx(q_hip_x)*([0; 0; l_upper_leg] + roty(q_knee_y)*[0;0; com_lower_leg]);
tip_position=rotz(q_hip_z)*roty(-q_hip_y)*rotx(q_hip_x)*([0; 0; l_upper_leg] + roty(q_knee_y)*[0;0; l_lower_leg]);

%transl(complex_leg.fkine([pi/6, pi/4, pi/2, -pi/2]))
%eval(subs(tip_position, [l_lower_leg, l_upper_leg, q_hip_z, q_hip_y, q_hip_x, q_knee_y], [L_lower_leg, L_upper_leg, pi/6, pi/4, pi/2, -pi/2]))

syms a1 a2 q1 q2 dq1 dq2

dr_upper_body = diff(r_upper_body, q_hip_z)*dq_hip_z + diff(r_upper_body, q_hip_y)*dq_hip_y + diff(r_upper_body, q_hip_x)*dq_hip_x; 
K_upper_body=1/2*mass_upper_leg*dr_upper_body'*dr_upper_body;
P_upper_body=mass_upper_leg*[0 0 g]*r_upper_body;

%simplify(eval(subs(K_upper_body, [l_lower_leg, l_upper_leg, q_hip_z, q_hip_y, q_hip_x, q_knee_y, dq_hip_z, dq_hip_y, dq_hip_x, dq_knee_y], [a1, a2, 0, q1, 0, q2, 0, dq1, 0, dq2])))
%simplify(eval(subs(P_upper_body, [l_lower_leg, l_upper_leg, q_hip_z, q_hip_y, q_hip_x, q_knee_y, dq_hip_z, dq_hip_y, dq_hip_x, dq_knee_y], [a1, a2, 0, q1, 0, q2, 0, dq1, 0, dq2])))

% lower body kinetic and potential energy
dr_lower_body = diff(r_lower_body, q_hip_z)*dq_hip_z + diff(r_lower_body, q_hip_y)*dq_hip_y + diff(r_lower_body, q_hip_x)*dq_hip_x + diff(r_lower_body, q_knee_y)*dq_knee_y; 
K_lower_body=1/2*mass_lower_leg*dr_lower_body'*dr_lower_body;
P_lower_body=mass_lower_leg*[0 0 g]*r_lower_body;

% upper body kinetic and potential energy
dr_upper_body = diff(r_upper_body, q_hip_z)*dq_hip_z + diff(r_upper_body, q_hip_y)*dq_hip_y + diff(r_upper_body, q_hip_x)*dq_hip_x; 
K_upper_body=1/2*mass_upper_leg*dr_upper_body'*dr_upper_body;
P_upper_body=mass_upper_leg*[0 0 g]*r_upper_body;

% langrangian
KE=K_upper_body + K_lower_body
PE=P_upper_body + P_lower_body;
L=K_upper_body + K_lower_body - P_upper_body - P_lower_body;

% gravity vector
G_vect = jacobian(PE,q).';
G_vect = simple(G_vect);

% mass-inertial matrix
D_mtx = simple(jacobian(KE,dq).');
D_mtx = simple(jacobian(D_mtx,dq));

% coriolis and centrifugal matrix
syms C_mtx real
n=max(size(q));
for k=1:n
	for j=1:n
		C_mtx(k,j)=0*g;
		for i=1:n
			C_mtx(k,j)=C_mtx(k,j)+1/2*(diff(D_mtx(k,j),q(i)) + ...
				diff(D_mtx(k,i),q(j)) - ...
				diff(D_mtx(i,j),q(k)))*dq(i);
		end
	end
end
C_mtx=simple(C_mtx);

%list_q = {}; 
%list_dq = {};
%param_list = {};
%write_fcn('D_matrix.m',{'q','params'},[list_q; param_list],{D_mtx,'D'});

phi=[ com_upper_leg, com_lower_leg, mass_upper_leg, mass_lower_leg, 1 ]';
W=simple(D_mtx*ddq+C_mtx*dq + G_vect);
%write_fcn('W_matrix.m',{'q','params'},[list_q; param_list],{W_mtx,'W'});

%subs(W,q,[0;0;0;0]);
%subs(subs(subs(W,q,[0;pi/2;0;0]),dq,[0;0;0;0]),ddq,[0;0;0;0])

com_adjustments=sym('ca', [5, 1]);
com_adjustments(5,1)=1;

mass_params=[mass_upper_leg; mass_lower_leg];
% 
% W_ca = zeros(4,5)
% 
% 
% [cw1, tw1]=coeffs(W(1),[com_upper_leg, com_lower_leg])
% % tw1 = [ com_upper_leg^2, com_lower_leg^2, com_lower_leg, 1]
% W_ca(1,1)=cw1(1,1);
% W_ca(1,2)= 0;
% W_ca(1,3)=cw1(1,2);
% W_ca(1,4)=cw1(1,3);
% W_ca(1,5)=cw1(1,4);
% 
% [cw2, tw2]=coeffs(W(2),[com_upper_leg, com_lower_leg])
% %tw2 = [ com_upper_leg^2, com_upper_leg, com_lower_leg^2, com_lower_leg, 1];
% W_ca(2,:)=cw2;
% 
% [cw3, tw3]=coeffs(W(3),[com_upper_leg, com_lower_leg])
% % tw3 = [ com_upper_leg^2, com_upper_leg, com_lower_leg^2, com_lower_leg, 1];
% W_ca(3,:)=cw3;
% 
% [cw4, tw4]=coeffs(W(4),[com_upper_leg, com_lower_leg])
% % tw4 = [ com_lower_leg^2, com_lower_leg];
% W_ca(1,3)=cw1(1,2);
% W_ca(1,4)=cw1(1,3);

% write_fcn('W_matrix.m',{'q','params'},[list_q; param_list],{W_mtx,'W'});

