
com_adjustments=sym('ca', [5, 1]);
com_adjustments(5,1)=1;
com_params=[com_upper_leg^2, com_upper_leg, com_lower_leg^2, com_lower_leg, 1];

mass_params=[mass_upper_leg; mass_lower_leg];

W_ca = sym('W_ca',[4,5])

[cw1, tw1]=coeffs(W(1),[com_upper_leg, com_lower_leg])
% tw1 = [ com_upper_leg^2, com_lower_leg^2, com_lower_leg, 1]
W_ca(1,1)=cw1(1,1);
W_ca(1,2)= 0;
W_ca(1,3)=cw1(1,2);
W_ca(1,4)=cw1(1,3);
W_ca(1,5)=cw1(1,4);

[cw2, tw2]=coeffs(W(2),[com_upper_leg, com_lower_leg])
%tw2 = [ com_upper_leg^2, com_upper_leg, com_lower_leg^2, com_lower_leg, 1];
W_ca(2,:)=cw2;

[cw3, tw3]=coeffs(W(3),[com_upper_leg, com_lower_leg])
% tw3 = [ com_upper_leg^2, com_upper_leg, com_lower_leg^2, com_lower_leg, 1];
W_ca(3,:)=cw3;

[cw4, tw4]=coeffs(W(4),[com_upper_leg, com_lower_leg])
% tw4 = [ com_lower_leg^2, com_lower_leg];
W_ca(4,3)=cw4(1,1);
W_ca(4,4)=cw4(1,2);

%(com_params'.*com_adjustments)
W_ca =  W_ca.*[com_params; com_params; com_params; com_params];

