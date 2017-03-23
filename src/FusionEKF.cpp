#include "FusionEKF.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;
    previous_timestamp_ = 0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
                0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
                0, 0.0009, 0,
                0, 0, 0.09;

    H_laser_<<1,0,0,0,
              0,1,0,0;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.P_ = MatrixXd(4, 4);

    VectorXd raw_measurements_ = measurement_pack.raw_measurements_;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = raw_measurements_(0);
      double phi = raw_measurements_(1);
      double rho_dot = raw_measurements_(2);
      ekf_.x_ << rho*cos(phi), rho*sin(phi), rho_dot*cos(phi), rho_dot*sin(phi);
      ekf_.P_ <<1,0,0,0,
                0,1,0,0,
                0,0,1,0,
                0,0,0,1;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << raw_measurements_(0), raw_measurements_(1), 0, 0;
      ekf_.P_ <<1,0,0,0,
                0,1,0,0,
                0,0,10,0,
                0,0,0,10;
    }

    previous_timestamp_ =  measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
   /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  //compute the time elapsed between the current and previous measurements
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_<< 1,0,dt,0,
            0,1,0,dt,
            0,0,1,0,
            0,0,0,1;

  double dt_2 = pow(dt,2);
  double dt_3 = pow(dt,3);
  double dt_4 = pow(dt,4);

  double noise_ax = 9;
  double noise_ay = 9;

  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_<< dt_4*noise_ax/4, 0, dt_3*noise_ax/2, 0,
            0, dt_4*noise_ay/4, 0, dt_3*noise_ay/2,
            dt_3*noise_ax/2, 0, dt_2*noise_ax, 0,
            0, dt_3*noise_ay/2, 0, dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {      
        // Radar updates
        Hj_ = tools.CalculateJacobian(ekf_.x_);
        if(!(Hj_.array()==1e9).all())
        {
            ekf_.H_ = Hj_;
            ekf_.R_ = R_radar_;
            ekf_.UpdateEKF(measurement_pack.raw_measurements_);
        }
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
