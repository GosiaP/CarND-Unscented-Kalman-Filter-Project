#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF()
  : std_a_(0.6)
  , std_yawdd_(0.3)
  , std_laspx_(0.15)
  , std_laspy_(0.15)
  , std_radr_(0.3)
  , std_radphi_(0.03)
  , std_radrd_(0.3)
  , n_x_(5)
  , n_aug_(7)
  , lambda_(3 - n_aug_)
  , n_z_radar_(3)
  , n_z_lidar_(2)
  , is_initialized_(false)
  , use_laser_(true)
  , use_radar_(true)
  , x_(VectorXd::Zero(n_x_))
  , x_aug_(n_aug_)
  , P_(MatrixXd::Identity(n_x_, n_x_))
  , P_aug_(n_aug_, n_aug_)
  , Xsig_aug_(n_aug_, 2 * n_aug_ + 1)
  , Xsig_pred_(n_x_, 2 * n_aug_ + 1)
  , previous_timestamp_(0)
  , NIS_radar_(0.)
  , NIS_lidar_(0.)
  , weights_(CalcWeights(lambda_, n_aug_))
  , R_radar_(CalcRadarMeasNoise(std_radr_, std_radphi_, std_radrd_))
  , R_lidar_(CalcLidarMeasNoise(std_laspx_, std_laspy_))
{
  std::cout
      << "R_radar_= " << std::endl << R_radar_ << std::endl
      << "R_lidar_= " << std::endl << R_lidar_ << std::endl;
}

UKF::~UKF() {}


VectorXd UKF::CalcWeights(double lambda, int n_aug)
{
  VectorXd w (2 * n_aug + 1);
  double weight_0 = lambda / (lambda + n_aug);
  w(0) = weight_0;

  //2n+1 weights
  for (int i = 1; i<2 * n_aug + 1; i++) 
  {  
    double weight = 0.5 / (n_aug + lambda);
    w(i) = weight;
  }
  return w;
}

MatrixXd UKF::CalcRadarMeasNoise(double std_radr, double std_radphi, double std_radrd)
{
  MatrixXd R (3, 3);
  R <<  std_radr * std_radr, 0., 0.,
        0., std_radphi*std_radphi, 0.,
        0., 0., std_radrd*std_radrd;

  return R;
}

MatrixXd UKF::CalcLidarMeasNoise(double std_laspx, double std_laspy)
{
  MatrixXd R (2, 2);
  R << std_laspx * std_laspx, 0.,
       0., std_laspy*std_laspy;

  return R;
}


/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) 
{
  /*****************************************************************************
  *  Initialization
  ****************************************************************************/
  if (!is_initialized_)
  {
    // first measurement
    previous_timestamp_ = meas_package.timestamp_;

    switch (meas_package.sensor_type_)
    {
      case MeasurementPackage::RADAR:
      {
        if (use_laser_)
        {
          /**
          Convert radar from polar to cartesian coordinates and initialize state.
          */
          double ro = meas_package.raw_measurements_(0);
          double phi = meas_package.raw_measurements_(1);
          double ro_dot = meas_package.raw_measurements_(2);

          x_(0) = ro * cos(phi);
          x_(1) = ro * sin(phi);
        }
      }
      break;

      case MeasurementPackage::LASER:
      {
        if (use_laser_)
        {
          /**
          Initialize state.
          */
          x_(0) = meas_package.raw_measurements_(0);
          x_(1) = meas_package.raw_measurements_(1);
        }
      }
      break;
    }

    previous_timestamp_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  
  /*****************************************************************************
  *  Prediction
  ****************************************************************************/
  //compute the time elapsed between the current and previous measurements
  auto dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;

  Prediction(dt);

  /*****************************************************************************
  *  Update
  ****************************************************************************/

  switch (meas_package.sensor_type_)
  {
    case MeasurementPackage::LASER:
    {
      UpdateLidar(meas_package);
    }
    break;

    case MeasurementPackage::RADAR:
    {
      UpdateRadar(meas_package);
    }
    break;
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) 
{
  /**
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  
  GenerateAugmentedSigmaPoints();
  PredictAugmentedSigmaPoints(delta_t);
  PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  //extract measurement as VectorXd
  const VectorXd& z = meas_package.raw_measurements_;
  VectorXd z_pred;
  MatrixXd Zsig, S;

  PredictLidar(z_pred, Zsig, S);
  UpdateLidar(z, z_pred, Zsig, S);
}

void UKF::PredictLidar(VectorXd& z_pred_out, MatrixXd& Zsig_out, MatrixXd& S_out)
{
  auto n_z = n_z_lidar_;
  MatrixXd Zsig{ n_z, 2 * n_aug_ + 1 };

  //transform sigma points into measurement space
  //2n+1 simga points
  for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  {  
    // extract values for better readibility
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);

    // measurement model
    Zsig(0, i) = p_x;
    Zsig(1, i) = p_y;
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  //2n+1 simga points
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {  
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix 
  S = S + R_lidar_;

  z_pred_out = z_pred;
  S_out = S;
  Zsig_out = Zsig;
}

void UKF::UpdateLidar(VectorXd const & z, VectorXd const & z_pred, MatrixXd const & Zsig, MatrixXd & S)
{
  //set measurement dimension, lidar can measure p_x and p_y
  const int n_z = n_z_lidar_;

  //create matrix for cross correlation Tc
  MatrixXd Tc{ n_x_, n_z };

  //calculate cross correlation matrix
  Tc.fill(0.0);

  //2n+1 simga points
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {  
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //calculate NIS
  NIS_lidar_ = z_diff.transpose() * S.inverse() * z_diff;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S*K.transpose();
}


void UKF::UpdateRadar(VectorXd const& z, VectorXd const& z_pred, MatrixXd const& Zsig, MatrixXd& S)
{
  const int n_z = n_z_radar_;

  //create matrix for cross correlation Tc
  MatrixXd Tc{ n_x_, n_z };

  //calculate cross correlation matrix
  Tc.fill(0.);

  //2n+1 simga points
  for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  {  
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = Tools::NormalizeAngle(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = Tools::NormalizeAngle(x_diff(3));
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;
  z_diff(1) = Tools::NormalizeAngle(z_diff(1));
  //calculate NIS
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S*K.transpose();
}

void UKF::PredictRadar(VectorXd & z_pred_out, MatrixXd & Zsig_out, MatrixXd & S_out)
{
  const int n_z = n_z_radar_;
  MatrixXd Zsig{ n_z, 2 * n_aug_ + 1 };

  //transform sigma points into measurement space
  //2n+1 simga points
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {  
    // extract values for better readibility
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v   = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0, i) = sqrt(p_x*p_x + p_y * p_y);                         //ro
    Zsig(1, i) = atan2(p_y, p_x);                                   //phi
    Zsig(2, i) = (p_x*v1 + p_y * v2) / sqrt(p_x*p_x + p_y * p_y);   //ro_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  //2n+1 simga points
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {  
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    z_diff(1) = Tools::NormalizeAngle(z_diff(1));
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S = S + R_radar_;

  z_pred_out = z_pred;
  S_out = S;
  Zsig_out = Zsig;
}


/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  //extract measurement as VectorXd
  VectorXd const& z = meas_package.raw_measurements_;
  VectorXd z_pred;
  MatrixXd Zsig, S;

  PredictRadar(z_pred, Zsig, S);
  UpdateRadar(z, z_pred, Zsig, S);
}

void UKF::GenerateAugmentedSigmaPoints()
{
  //create augmented mean state
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  //create augmented covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5, 5) = P_;
  P_aug_(5, 5) = std_a_ * std_a_;
  P_aug_(6, 6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug_.llt().matrixL();

  //create augmented sigma points
  Xsig_aug_.col(0) = x_aug_;
  for (int i = 0; i < n_aug_; i++) 
  {
    Xsig_aug_.col(i + 1) = x_aug_ + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug_.col(i + 1 + n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L.col(i);
  }
}

void UKF::PredictAugmentedSigmaPoints(double delta_t)
{
  //predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    //extract values for better readability
    double p_x      = Xsig_aug_(0, i);
    double p_y      = Xsig_aug_(1, i);
    double v        = Xsig_aug_(2, i);
    double yaw      = Xsig_aug_(3, i);
    double yawd     = Xsig_aug_(4, i);
    double nu_a     = Xsig_aug_(5, i);
    double nu_yawdd = Xsig_aug_(6, i);
      
    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001)
    {
      px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else 
    {
      px_p = p_x + v * delta_t*cos(yaw);
      py_p = p_y + v * delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p  = v_p + nu_a * delta_t;

    yaw_p  = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }
}

void UKF::PredictMeanAndCovariance()
{
  //predicted state mean
  x_.fill(0.0);
  //iterate over sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {  
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  //iterate over sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {  
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = Tools::NormalizeAngle(x_diff(3));
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

