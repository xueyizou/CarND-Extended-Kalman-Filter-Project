#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

const double PI = 3.1415926535897932384626433832795;
void normalize_phi(double & phi)
{
    while(phi<-PI)
    {
        phi += 2*PI;
    }

    while(phi>PI)
    {
        phi -= 2*PI;
    }
}

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict() {
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
    normalize_phi(y(1));
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

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    double rho = sqrt(x_(0)*x_(0)+x_(1)*x_(1));
    double phi = atan2(x_(1), x_(0));
    double rho_dot = (x_(2)*x_(0) + x_(3)*x_(1))/rho;

    VectorXd z_pred(3);
    z_pred<<rho, phi, rho_dot;

    VectorXd y = z - z_pred;
    normalize_phi(y(1));
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
