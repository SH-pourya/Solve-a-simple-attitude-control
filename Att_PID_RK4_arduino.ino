const double Tstep = 0.0001;    // Timestep used by the algorithm
const int n = 12;              // number of state variables in the ODE
//int steps=0;                  // number of calculated steps
double t = 0;               // current time
double pi = 3.1415926536;
double time_final = 20;
/*
 
  X[6] = omega_x ,  X[7] = omega_y , X[8] = omega_z , X[9] = phi , X[10] = theta , X[11] = psi ,
*/
double X[n] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
double X_new[n] = {};
double Y[n] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
double XDOT[n] = {};

double K1[n], K2[n], K3[n], K4[n];   // arrays for the Runge-Kutta intermediate points

void setup() {
  // put your setup code here, to run once:
  Serial.begin(57600);
}

void loop() {
  Serial.print("time: ");
  Serial.print(t);
  Serial.print(",\t");
  Serial.print("phi: ");
  Serial.print(X[9]*180/pi,5);
  Serial.print(",\t");
  Serial.print("theta: ");
  Serial.print(X[10]*180/pi,5);
  Serial.print(",\t");
  Serial.print("psi: ");
  Serial.println(X[11]*180/pi,5);
//  Serial.println("------------------");

  if (t >= time_final) {
    for (;;);                  //STOP LOOP
  }


  for (int i = 0; i < n; i++) {
    Y[i] = X[i];
  }
  my_function(Y, t , XDOT);
  for (int i = 0; i < n; i++) {
    K1[i] = Tstep * XDOT[i];
    Y[i] = X[i] + 0.5 * K1[i];
  }

  my_function(Y, t + Tstep / 2 , XDOT);
  for (int i = 0; i < n; i++) {
    K2[i] = Tstep * XDOT[i];
    Y[i] = X[i] + 0.5 * K2[i];
  }

  my_function(Y, t + Tstep / 2 , XDOT);
  for (int i = 0; i < n; i++) {
    K3[i] = Tstep * XDOT[i];
    Y[i] = X[i] + K3[i];
  }

  my_function(Y, t + Tstep , XDOT);
  for (int i = 0; i < n; i++) {
    K4[i] = Tstep * XDOT[i];
  }

  for (int i = 0; i < n; i++) {
    //    X_new[i] = X[i] + ((1 / 6) * (K1[i] + (2 * K2[i]) + (2 * K3[i]) + K4[i]));
    X_new[i] = X[i] +  (K1[i] + (2 * K2[i]) + (2 * K3[i]) + K4[i]) / 6;
  }

  for (int i = 0; i < n; i++) {
    X[i] = X_new[i];
  }

  t = t + Tstep;
  //  Serial.println(X[0]);
  //  Serial.println(X[1]);
  //  Serial.println(X[2]);
  //  Serial.println("------------------");

  //  delay(1000);
}


void my_function(double Y[12] , double t , double outV[12]) {
  double I[3][3] = {{200, 0, 0}, {0, 150, 0}, {0, 0, 300}};
  double att_ref[] = {0.5236, 0.5236, 0.5236};
  double PID_Xi[] = {220.3, 27.536, 376, 100}; 
  double PID_Yi[] = {165.2, 20.652, 282, 100}; 
  double PID_Zi[] = {330.4, 41.304, 564, 100};
//  double R_e = 6378000;
//  double omega = 0.000072921;
  double g = 9.81;
  double Ixx = I[0][0];
  double Iyy = I[1][1];
  double Izz = I[2][2];
//  double Ixxy = I[0][1];
//  double Ixxz = I[0][2];
//  double Iyyx = I[1][0];
//  double Iyyz = I[1][2];
//  double Izzx = I[2][0];
//  double Izzy = I[2][1];
  double phi_refrence = att_ref[0];
  double theta_refrence = att_ref[1];
  double psi_refrence = att_ref[2];
    
  double Kp_x = PID_Xi[0];
  double Ki_x = PID_Xi[1];
  double Kd_x = PID_Xi[2];
  double N_x = PID_Xi[3];
    
  double Kp_y = PID_Yi[0];
  double Ki_y = PID_Yi[1];
  double Kd_y = PID_Yi[2]; 
  double N_y = PID_Yi[3];
    
  double Kp_z = PID_Zi[0];
  double Ki_z = PID_Zi[1];
  double Kd_z = PID_Zi[2];    
  double N_z = PID_Zi[3];
  
  double x2_phi = Y[0]; 
  double xphi_feed = Y[1]; 
  double x2_theta = Y[2]; 
  double xtheta_feed = Y[3]; 
  double x2_psi = Y[4]; 
  double xpsi_feed = Y[5]; 
  double omega_x = Y[6]; 
  double omega_y = Y[7]; 
  double omega_z = Y[8]; 
  double phi = Y[9]; 
  double theta = Y[10]; 
  double psi = Y[11];

    //--------------- calculate of Mx ----------------------
    double error_phi = phi_refrence - phi;
    double x1_phi = Kp_x * error_phi;
    double x3_phi = (Kd_x * (error_phi) - xphi_feed) * N_x;
    double x2_phiD = Ki_x * error_phi;
    double xphi_feedD = x3_phi;
    double Mx = x1_phi + x2_phi + x3_phi;
    //--------------- calculate of My ----------------------
    double error_theta = theta_refrence - theta;
    double x1_theta = Kp_y * error_theta;
    double x3_theta = (Kd_y * (error_theta) - xtheta_feed) * N_y;
    double x2_thetaD = Ki_y * error_theta;
    double xtheta_feedD = x3_theta;
    double My = x1_theta + x2_theta + x3_theta;
    //--------------- calculate of Mz ----------------------
    double error_psi = psi_refrence - psi;
    double x1_psi = Kp_z * error_psi;
    double x3_psi = (Kd_z * (error_psi) - xpsi_feed) * N_z;
    double x2_psiD = Ki_z * error_psi;
    double xpsi_feedD = x3_psi;
    double Mz = x1_psi + x2_psi + x3_psi;    

    double omegaD_x = (Mx - (omega_y * omega_z * (Izz - Iyy))) / Ixx;
    double omegaD_y = (My - (omega_x * omega_z * (Ixx - Izz))) / Iyy;
    double omegaD_z = (Mz - (omega_x * omega_y * (Iyy - Ixx))) / Izz;

    double p = omega_x; 
    double q = omega_y; 
    double r = omega_z;
    double phiD = p + ((q * sin(phi) + r * cos(phi)) * tan(theta));
    double thetaD = q * cos(phi) - r * sin(phi);
    double psiD = (q * sin(phi) + r * cos(phi)) * (1 / cos(theta));

    outV[0] = x2_phiD; outV[1] = xphi_feedD; outV[2] = x2_thetaD; 
    outV[3] = xtheta_feedD; outV[4] = x2_psiD; outV[5] = xpsi_feedD; 
    outV[6] = omegaD_x; outV[7] =omegaD_y; outV[8] = omegaD_z; 
    outV[9] = phiD; outV[10] = thetaD; outV[11] = psiD;
}
