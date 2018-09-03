#ifndef __FEM_INTERFACE_H__
#define __FEM_INTERFACE_H__
#include "../fem/fem.h"
#include "GLfunctions.h"

namespace GUI
{
void DrawWorld(const std::shared_ptr<FEM::World>& world);
void DrawConstraint(const std::shared_ptr<FEM::Cst>& c,const Eigen::VectorXd& x);

};

#endif
