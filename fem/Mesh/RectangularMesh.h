#ifndef __RECTANGULAR_MESH_H__
#define __RECTANGULAR_MESH_H__

#include "Mesh.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
namespace FEM
{

class RectangularMesh : public Mesh
{
public:
	RectangularMesh(const RectangularMesh& other) = delete;
	RectangularMesh& operator=(const RectangularMesh& other) = delete;
	std::shared_ptr<Mesh> Clone() override;
	static std::shared_ptr<RectangularMesh> Create(double _x = 1.0,double _y = 1.0, double _z = 1.0, int _nx = 5,int _ny = 5, int _nz = 5,const Eigen::Affine3d& M = Eigen::Affine3d::Identity());
protected:
	RectangularMesh(double _x = 1.0,double _y = 1.0, double _z = 1.0, int _nx = 5,int _ny = 5, int _nz = 5,const Eigen::Affine3d& _M = Eigen::Affine3d::Identity());

	int mNx,mNy,mNz;
	double mX,mY,mZ;
	Eigen::Affine3d mT;
};
};
#endif