#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF 
{
public:

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  const double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  const double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  const double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  const double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  const double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  const double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  const double std_radrd_;

  ///* State dimension
  const int n_x_;

  ///* Augmented state dimension
  const int n_aug_;

  ///* Sigma point spreading parameter
  const double lambda_;

  ///* Radar measurement dimension
  const int n_z_radar_;

  ///* Lidar measurement dimension
  const int n_z_lidar_;

  ///** initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///** if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* augmented state vector
  VectorXd x_aug_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* augmented state covariance matrix
  MatrixXd P_aug_;

  ///* predicted augmented sigma points matrix
  MatrixXd Xsig_aug_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long previous_timestamp_;

  ///* Normalized Innovation Squared (NIS) for radar measurements
  double NIS_radar_;

  ///* Normalized Innovation Squared (NIS) for lidar measurements
  double NIS_lidar_;

  ///* Vector for weights of sigma points
  const VectorXd weights_;

  ///* Measurement noise matrix radar
  const MatrixXd R_radar_;

  ///* Measurement noise matrix laser
  const MatrixXd R_lidar_;


  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

private:
  static VectorXd CalcWeights(double lambda, int n_aug);
  static MatrixXd CalcRadarMeasNoise(double std_radr, double std_radphi, double std_radrd);
  static MatrixXd CalcLidarMeasNoise(double std_laspx, double std_laspy);

  void GenerateAugmentedSigmaPoints();
  void PredictAugmentedSigmaPoints(double delta_t);
  void PredictMeanAndCovariance();

  /* Helper methods */
  void PredictRadar(VectorXd& z_pred_out, MatrixXd& Zsig_out, MatrixXd& S_out);
  void UpdateRadar(VectorXd const& z, VectorXd const & z_pred, MatrixXd const & Zsig, MatrixXd & S);
  void PredictLidar(VectorXd& z_pred_out, MatrixXd& Zsig_out, MatrixXd& S_out);
  void UpdateLidar(VectorXd const& z, VectorXd const & z_pred, MatrixXd const & Zsig, MatrixXd & S);
};

#endif /* UKF_H */
