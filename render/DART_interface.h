#ifndef __DART_INTERFACE_H__
#define __DART_INTERFACE_H__
#include "dart/dart.hpp"
#include "dart/gui/gui.hpp"
#include "dart/math/math.hpp"
#include "dart/simulation/simulation.hpp"
#include "GLfunctions.h"
typedef std::pair<dart::dynamics::BodyNode*,Eigen::Vector3d> AnchorPoint;
namespace GUI
{
void DrawSkeleton(
	const dart::dynamics::SkeletonPtr& skel,
	const Eigen::Vector3d& color=Eigen::Vector3d(0.8,0.8,0.8),bool box = true);

void DrawShape(const Eigen::Isometry3d& T,
	const dart::dynamics::Shape* shape,
	const Eigen::Vector3d& color);

void DrawMuscleWayPoints(const std::vector<AnchorPoint>& ap);
};


#endif