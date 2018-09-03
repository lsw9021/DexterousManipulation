#ifndef __MESH_H__
#define __MESH_H__
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <memory>
namespace FEM
{
class Mesh
{
public:
	~Mesh(){mVertices.clear();mTetrahedrons.clear();};
	virtual const std::vector<Eigen::Vector3d>& GetVertices(){return mVertices;};
	virtual const std::vector<Eigen::Vector4i>& GetTetrahedrons(){return mTetrahedrons;};

	Mesh(const Mesh& other) = delete;
	Mesh& operator=(const Mesh& other) = delete;
	virtual std::shared_ptr<Mesh> Clone() = 0;
protected:
	Mesh(){};
	std::vector<Eigen::Vector3d> mVertices;
	std::vector<Eigen::Vector4i> mTetrahedrons;
};
};


#endif