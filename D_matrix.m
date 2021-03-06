function [D] = D_matrix(q,params)

  D(1,1)=com_upper_leg^2*mass_upper_leg + l_upper_leg^2*mass_lower_leg + com_lower_leg^2*...
         mass_lower_leg*cos(q_hip_y)^2 + com_lower_leg^2*mass_lower_leg*cos(q_knee_y)^2 - l_upper_leg^2*mass_lower_leg*...
         cos(q_hip_x)^2*cos(q_hip_y)^2 + 2*com_lower_leg*l_upper_leg*mass_lower_leg*cos(q_knee_y) - com_upper_leg^2*...
         mass_upper_leg*cos(q_hip_x)^2*cos(q_hip_y)^2 - com_lower_leg^2*mass_lower_leg*cos(q_hip_y)^2*cos(q_knee_y)^2 -...
          com_lower_leg^2*mass_lower_leg*cos(q_hip_x)^2*cos(q_hip_y)^2*cos(q_knee_y)^2 - 2*com_lower_leg*l_upper_leg*...
         mass_lower_leg*cos(q_hip_x)^2*cos(q_hip_y)^2*cos(q_knee_y) - 2*com_lower_leg*l_upper_leg*mass_lower_leg*...
         cos(q_hip_x)*cos(q_hip_y)*sin(q_hip_y)*sin(q_knee_y) - 2*com_lower_leg^2*mass_lower_leg*cos(q_hip_x)*...
         cos(q_hip_y)*cos(q_knee_y)*sin(q_hip_y)*sin(q_knee_y);
  D(1,2)=-sin(q_hip_x)*(com_upper_leg^2*mass_upper_leg*cos(q_hip_x)*cos(q_hip_y) + l_upper_leg^2*...
         mass_lower_leg*cos(q_hip_x)*cos(q_hip_y) + com_lower_leg^2*mass_lower_leg*cos(q_knee_y)*sin(q_hip_y)*...
         sin(q_knee_y) + com_lower_leg*l_upper_leg*mass_lower_leg*sin(q_hip_y)*sin(q_knee_y) + com_lower_leg^2*...
         mass_lower_leg*cos(q_hip_x)*cos(q_hip_y)*cos(q_knee_y)^2 + 2*com_lower_leg*l_upper_leg*mass_lower_leg*cos(q_hip_x)*...
         cos(q_hip_y)*cos(q_knee_y));
  D(1,3)=com_upper_leg^2*mass_upper_leg*sin(q_hip_y) + l_upper_leg^2*mass_lower_leg*sin(q_hip_y) +...
          com_lower_leg^2*mass_lower_leg*cos(q_knee_y)^2*sin(q_hip_y) + 2*com_lower_leg*l_upper_leg*mass_lower_leg*...
         cos(q_knee_y)*sin(q_hip_y) - com_lower_leg*l_upper_leg*mass_lower_leg*cos(q_hip_x)*cos(q_hip_y)*sin(q_knee_y) -...
          com_lower_leg^2*mass_lower_leg*cos(q_hip_x)*cos(q_hip_y)*cos(q_knee_y)*sin(q_knee_y);
  D(1,4)=com_lower_leg*mass_lower_leg*cos(q_hip_y)*sin(q_hip_x)*(com_lower_leg + l_upper_leg*cos(q_knee_y));
  D(2,1)=-sin(q_hip_x)*(com_upper_leg^2*mass_upper_leg*cos(q_hip_x)*cos(q_hip_y) + l_upper_leg^2*...
         mass_lower_leg*cos(q_hip_x)*cos(q_hip_y) + com_lower_leg^2*mass_lower_leg*cos(q_knee_y)*sin(q_hip_y)*...
         sin(q_knee_y) + com_lower_leg*l_upper_leg*mass_lower_leg*sin(q_hip_y)*sin(q_knee_y) + com_lower_leg^2*...
         mass_lower_leg*cos(q_hip_x)*cos(q_hip_y)*cos(q_knee_y)^2 + 2*com_lower_leg*l_upper_leg*mass_lower_leg*cos(q_hip_x)*...
         cos(q_hip_y)*cos(q_knee_y));
  D(2,2)=com_lower_leg^2*mass_lower_leg + com_upper_leg^2*mass_upper_leg*cos(q_hip_x)^2 -...
          com_lower_leg^2*mass_lower_leg*cos(q_knee_y)^2 + l_upper_leg^2*mass_lower_leg*cos(q_hip_x)^2 + com_lower_leg^2*...
         mass_lower_leg*cos(q_hip_x)^2*cos(q_knee_y)^2 + 2*com_lower_leg*l_upper_leg*mass_lower_leg*cos(q_hip_x)^2*cos(q_knee_y);
  D(2,3)=-com_lower_leg*mass_lower_leg*sin(q_hip_x)*sin(q_knee_y)*(l_upper_leg + com_lower_leg*cos(q_knee_y));
  D(2,4)=-com_lower_leg*mass_lower_leg*cos(q_hip_x)*(com_lower_leg + l_upper_leg*cos(q_knee_y));
  D(3,1)=com_upper_leg^2*mass_upper_leg*sin(q_hip_y) + l_upper_leg^2*mass_lower_leg*sin(q_hip_y) +...
          com_lower_leg^2*mass_lower_leg*cos(q_knee_y)^2*sin(q_hip_y) + 2*com_lower_leg*l_upper_leg*mass_lower_leg*...
         cos(q_knee_y)*sin(q_hip_y) - com_lower_leg*l_upper_leg*mass_lower_leg*cos(q_hip_x)*cos(q_hip_y)*sin(q_knee_y) -...
          com_lower_leg^2*mass_lower_leg*cos(q_hip_x)*cos(q_hip_y)*cos(q_knee_y)*sin(q_knee_y);
  D(3,2)=-com_lower_leg*mass_lower_leg*sin(q_hip_x)*sin(q_knee_y)*(l_upper_leg + com_lower_leg*cos(q_knee_y));
  D(3,3)=com_upper_leg^2*mass_upper_leg + l_upper_leg^2*mass_lower_leg + com_lower_leg^2*...
         mass_lower_leg*cos(q_knee_y)^2 + 2*com_lower_leg*l_upper_leg*mass_lower_leg*cos(q_knee_y);
  D(3,4)=0;
  D(4,1)=com_lower_leg*mass_lower_leg*cos(q_hip_y)*sin(q_hip_x)*(com_lower_leg + l_upper_leg*cos(q_knee_y));
  D(4,2)=-com_lower_leg*mass_lower_leg*cos(q_hip_x)*(com_lower_leg + l_upper_leg*cos(q_knee_y));
  D(4,3)=0;
  D(4,4)=com_lower_leg^2*mass_lower_leg;

 