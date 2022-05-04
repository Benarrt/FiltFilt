#ifndef EIGEN_UTIL_H
#define EIGEN_UTIL_H

#include <vector>
#include "eigen3/Eigen/Dense"

namespace Util
{
    Eigen::VectorXd toEigen(const std::vector<double>& other, size_t fromElement = 0, size_t toElement = 0)
    {
        if(toElement == 0 || toElement > other.size())
            toElement = other.size();

        if(fromElement > other.size())
            fromElement = other.size();

        Eigen::VectorXd vec(toElement-fromElement);
        for(int i = fromElement; i < toElement; ++i)
            vec(i - fromElement) = other[i];

        return vec;
    }
};

#endif //EIGEN_UTIL_H