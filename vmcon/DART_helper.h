#ifndef __DART_HELPER_H__
#define __DART_HELPER_H__
#include "dart/dart.hpp"
#include "dart/gui/gui.hpp"
#include "dart/math/math.hpp"
#include "dart/simulation/simulation.hpp"

enum JOINT_TYPE
{
	WELD,
	REVOLUTE,
	UNIVERSAL,
	EULER,
	BALL_AND_SOCKET,
	FREE
};
void MakeRootBody(
	const dart::dynamics::SkeletonPtr& skel,
	const std::string& name,
	const Eigen::Vector3d& size,
	const Eigen::Vector3d& c_to_joint,
	JOINT_TYPE joint_type,
	double mass);

void MakeBody(
	const dart::dynamics::SkeletonPtr& skel,
	const dart::dynamics::BodyNodePtr& parent,
	const std::string& name,
	const Eigen::Vector3d& size,
	const Eigen::Vector3d& p_to_joint,
	const Eigen::Vector3d& c_to_joint,
	JOINT_TYPE joint_type,
	double mass);

void MakeBall(
	const dart::dynamics::SkeletonPtr& skel,
	const Eigen::Vector3d& init_pos,
	double rad,
	double mass);
void MakeDumbell(
		const dart::dynamics::SkeletonPtr& skel,
	const Eigen::Vector3d& init_pos,
	double rad,
	double mass);
#endif