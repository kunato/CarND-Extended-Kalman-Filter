#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
   /**
   * TODO: Calculate the RMSE here.
   */
   VectorXd result = VectorXd(4);
   result << 0, 0, 0, 0;

   for (unsigned int i = 0; i < estimations.size(); ++i)
   {

      VectorXd residual = estimations[i] - ground_truth[i];

      // coefficient-wise multiplication
      residual = residual.array() * residual.array();
      result += residual;
   }
   result = result / estimations.size();

   result = result.array().sqrt();

   return result;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state)
{
   /**
   * Done:
   * Calculate a Jacobian here.
   */
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // pre-compute a set of terms to avoid repeated calculation
  float c1 = px * px + py * py;
  float c2 = sqrt(c1);
  float c3 = (c1 * c2);

  MatrixXd Hj = MatrixXd(3, 4);
  // check division by zero
  if (fabs(c1) < 0.0001)
  {
    return Hj;
  }
  else
  {
    // compute the Jacobian matrix
    Hj << (px / c2), (py / c2), 0, 0,
        -(py / c1), (px / c1), 0, 0,
        py * (vx * py - vy * px) / c3, px * (px * vy - py * vx) / c3, px / c2, py / c2;
        return Hj;
  }
}

VectorXd Tools::Cart2Polar(const VectorXd& x_state)
{
	//recover state parameters
	double px = x_state(0);
	double py = x_state(1);
	double vx = x_state(2);
	double vy = x_state(3);

	double rho = sqrt(px*px + py*py);
	double phi = atan2(py, px);
	double rho_dot = (rho != 0) ? (px*vx +py*vy) / rho : 0;

	VectorXd x_polar(3);
	x_polar << rho, phi, rho_dot;

	return x_polar;
}