#ifndef __OBJ_MESH_H__
#define __OBJ_MESH_H__

#include "Mesh.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>

namespace FEM
{

class OBJMesh : public Mesh
{
public:
	OBJMesh(const OBJMesh& other) = delete;
	OBJMesh& operator=(const OBJMesh& other) = delete;
	std::shared_ptr<Mesh> Clone() override;
	static std::shared_ptr<OBJMesh> Create(const std::string& path,const Eigen::Affine3d& M = Eigen::Affine3d::Identity());
protected:
	OBJMesh(const std::string& path,const Eigen::Affine3d& _M = Eigen::Affine3d::Identity());

	std::string mPath;
	Eigen::Affine3d mT;
};
};
#endif