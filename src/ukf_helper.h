#ifndef UKFHelper_H
#define UKFHelper_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKFHelper {
public:

  /**
   * Constructor
   */
  UKFHelper();

  /**
   * Destructor
   */
  virtual ~UKFHelper();

};

#endif /* UKFHelper_H */
