//
// Created by 汪子琦 on 3/6/17.
//

#ifndef SUPPORTER_NORMALIZING_MODEL_H
#define SUPPORTER_NORMALIZING_MODEL_H

#include "Eigen/Core"
#include "settings.h"

class NormalizingModel
{
public:
    NormalizingModel()
    {

    }
public:

    void size_normalize(Eigen::MatrixXd &V)
    {
        Eigen::Vector3d m = V.colwise().minCoeff();
        Eigen::Vector3d M = V.colwise().maxCoeff();

        //scale
        //V  *= (settings_.maximum_model_heigh / (M(1) - m(1)));

        m = V.colwise().minCoeff();
        M = V.colwise().maxCoeff();
        for(size_t kd = 0; kd < V.rows(); kd++)
        {

            //transform
            V.row(kd) -=  Eigen::RowVector3d((m(0) + M(0)) / 2 , m(1), (m(2) + M(2)) / 2 );
            V.row(kd) +=  Eigen::RowVector3d(11.0 / 2 * settings_.pad_size, 0, 9.0 / 2  * settings_.pad_size);

        }
    }

private:

    Settings settings_;
};

#endif //SUPPORTER_NORMALIZING_MODEL_H
