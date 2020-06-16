
#include "ukf.h"
#include "Eigen/Dense"
#include<iostream>
#include<stdio.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // state vector dimension
  n_x_ = 5;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 0.0225, 0,
        0, 0, 0, 0, 0.0225;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.8;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
   * End DO NOT MODIFY section for measurement noise values
   */

  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */

   // Augmented state dimension
   n_aug_ = n_x_ + 2;

   lambda_ = 3 - n_aug_;

   time_us_ = 0.0;

   // Initialize weights_
   weights_ = VectorXd(2 * n_aug_ + 1);
   weights_.fill(0.5/(lambda_ + n_aug_));
   weights_(0) = lambda_/(lambda_ + n_aug_);

   // sigma point prediction
   Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

   // Initialize measurement noise covairance matix
   R_radar_ = MatrixXd(3, 3);
   R_radar_ << std_radr_*std_radr_, 0, 0,
               0, std_radphi_*std_radphi_, 0,
               0, 0, std_radrd_*std_radrd_;

   R_lidar_ = MatrixXd(2,2);
   R_lidar_ << std_laspx_*std_laspx_, 0,
               0, std_laspy_*std_laspy_;

   // the current NIS for radar
   NIS_radar_ = 0.0;

   // the current NIS for laser
   NIS_laser_ = 0.0;
}

UKF::~UKF() {}
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
 if(! is_initialized_){
   if(meas_package.sensor_type_==MeasurementPackage::LASER){
    x_<<meas_package.raw_measurements_[0],meas_package.raw_measurements_[1],0,0,0;
    
  }
  else if (meas_package.sensor_type_==MeasurementPackage::RADAR){
    double range=meas_package.raw_measurements_[0];
    double phi=meas_package.raw_measurements_[1];
    double RR=meas_package.raw_measurements_[2];
    double x = range * cos(phi);
    double y = range * sin(phi);
    double vx = RR * cos(phi);
    double vy = RR * sin(phi);
    double v = sqrt(vx * vx + vy * vy);
    x_<<x,y,v,0,0;
  }
   is_initialized_=true;
   return;
 }
   double dt = (meas_package.timestamp_ -time_us_)/1000000.0;
   time_us_ = meas_package.timestamp_;
   Prediction(dt);
 if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
     UpdateRadar(meas_package);
   }

   if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
     UpdateLidar(meas_package);
   }
}
void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  /* initialize augmented sigman points*/
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  x_aug.head(5)=x_;
  x_aug(5)=0;
  x_aug(6)=0;
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;
  MatrixXd L = P_aug.llt().matrixL();
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; ++i) {
    Xsig_aug.col(i+1)= x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  /*predict sigma points*/
  
  for(int j=0;j<2*n_aug_+1;j++){
    double p_x=Xsig_aug(0,j);
    double p_y=Xsig_aug(1,j);
    double v = Xsig_aug(2,j);
    double yaw = Xsig_aug(3,j);
    double yawd = Xsig_aug(4,j);
    double nu_a = Xsig_aug(5,j);
    double nu_yawdd = Xsig_aug(6,j);
    double px_p, py_p;

    // avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * (-cos(yaw + yawd*delta_t) + cos(yaw));
    } else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    // write predicted sigma point into right column
    Xsig_pred_(0,j) = px_p;
    Xsig_pred_(1,j) = py_p;
    Xsig_pred_(2,j) = v_p;
    Xsig_pred_(3,j) = yaw_p;
    Xsig_pred_(4,j) = yawd_p;
    
  }
  
  
  x_.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++){
   x_=x_+weights_(i) * Xsig_pred_.col(i); 
  }
  P_.fill(0.0);
   for(int i = 0; i < 2*n_aug_+1; i++){
     VectorXd x_diff = Xsig_pred_.col(i) - x_;
     // angle normalization
     while(x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
     while(x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;
     P_ = P_ + weights_(i)*x_diff*x_diff.transpose();
   }
  
  
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  VectorXd z=meas_package.raw_measurements_;
  
  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(2, 2 * n_aug_ + 1);
  // mean predicted measurement
  VectorXd z_pred = VectorXd(2);
  // measurement covariance matrix S
  MatrixXd S = MatrixXd(2,2);
  
  for(int i=0;i<2*n_aug_+1;i++){
    Zsig(0,i)=Xsig_pred_(0,i);
    Zsig(1,i)=Xsig_pred_(1,i);
  }
  
  z_pred.fill(0.0);
  //predict mean
  for (int i=0; i < 2*n_aug_+1; ++i) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  //predict coverience
  S.fill(0.0);
  for(int i = 0; i < 2*n_aug_+1; i++){
     VectorXd z_diff = Zsig.col(i) - z_pred;

     S = S + weights_(i) * z_diff * z_diff.transpose();
   }
  // add measurement noise covariance matrix
  
 
 
   S = S + R_lidar_;
  //calculate cross correlation
  MatrixXd Tc = MatrixXd(n_x_, 2);
   Tc.fill(0.0);
   for(int i = 0; i < 2*n_aug_+1; i++){
     VectorXd x_diff = Xsig_pred_.col(i) - x_;

     VectorXd z_diff = Zsig.col(i) - z_pred;

     Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
   }
  // calculate Kalman gain K
  MatrixXd K = Tc * S.inverse();
  VectorXd z_diff = z - z_pred;
   x_ = x_ + K * z_diff;
   P_ = P_ - K*S*K.transpose();
  
  //calculate NIS
   NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
  
  
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  VectorXd z_ = meas_package.raw_measurements_;
  MatrixXd Zsig = MatrixXd(3, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(3);
  MatrixXd S = MatrixXd(3,3);
  
  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // extract values for better readability
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double yawd = Xsig_pred_(4, i);
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                       // r
    Zsig(1,i) = atan2(p_y,p_x);                                // phi
    Zsig(2,i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   // r_dot
  }
  S.fill(0.0);
  z_pred.fill(0.0);
   for(int i = 0; i < 2*n_aug_+1; i++){
     z_pred = z_pred + weights_(i)*Zsig.col(i);
   }
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  
  
  S = S + R_radar_;
  MatrixXd Tc = MatrixXd(n_x_, 3);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
MatrixXd K = Tc * S.inverse();
  
  VectorXd z_diff = z_ - z_pred;
  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  //calculate NIS
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
  
}