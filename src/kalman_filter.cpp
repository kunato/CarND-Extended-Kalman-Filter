#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}
// use F, Q || P, x
// State transition, un-certentainty cov state, noise matrix
void KalmanFilter::Predict()
{
  /**
   * DONE: predict the state
   */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

// use H, R || P, x
// un-certentainty cov state, measurement->state func, measurement noise
void KalmanFilter::Update(const VectorXd &z)
{
  // z base on measurement func.
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{
  /**
   * update the state by using Extended Kalman Filter equations
   */

  MatrixXd Hj = tools.CalculateJacobian(x_);

  // z base on polar coord
  // we'll use the h function directly to map predicted locations x'x from Cartesian to polar coordinates.
  VectorXd z_pred = tools.Cart2Polar(x_);
  VectorXd y = z - z_pred;

	double phi = y(1);
	if (phi > M_PI){
    y(1) = phi - 2 * M_PI;
  }
	else if (phi < -M_PI){
    y(1) = phi + 2 * M_PI;
  }

  // after this will same as normal Update with Hj is used instead of H.
  MatrixXd Ht = Hj.transpose();
  MatrixXd S = Hj * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj) * P_;
}
