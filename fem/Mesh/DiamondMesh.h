#ifndef __DIAMOND_MESH_H__
#define __DIAMOND_MESH_H__

#include "RectangularMesh.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
namespace FEM
{

class DiamondMesh : public RectangularMesh
{
public:
	int GetStartingPointIndex(){return mStartPointIndex;};
	int GetEndingPointIndex(){return mEndPointIndex;};

	DiamondMesh(const DiamondMesh& other) = delete;
	DiamondMesh& operator=(const DiamondMesh& other) = delete;
	std::shared_ptr<Mesh> Clone() override;
	static std::shared_ptr<DiamondMesh> Create(double _x = 1.0,double _y = 1.0, double _z = 1.0, int _nx = 5,int _ny = 5, int _nz = 5,const Eigen::Affine3d& M = Eigen::Affine3d::Identity());
protected:
	DiamondMesh(double _x = 1.0,double _y = 1.0, double _z = 1.0, int _nx = 5,int _ny = 5, int _nz = 5,const Eigen::Affine3d& M = Eigen::Affine3d::Identity());

	int mStartPointIndex,mEndPointIndex;


};
};

#endif