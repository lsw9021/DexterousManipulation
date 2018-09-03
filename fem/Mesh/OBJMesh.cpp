#include "OBJMesh.h"
#include <fstream>
#include <sstream>
#include <iostream>
using namespace FEM;
std::shared_ptr<Mesh>
OBJMesh::
Clone()
{
	return OBJMesh::Create(mPath,mT);
}
std::shared_ptr<OBJMesh>
OBJMesh::
Create(const std::string& path,const Eigen::Affine3d& _M)
{
	auto rm = new OBJMesh(path,_M);

	return std::shared_ptr<OBJMesh>(rm);
}
OBJMesh::
OBJMesh(const std::string& path,const Eigen::Affine3d& _M)
	:mPath(path),mT(_M)
{
	std::ifstream ifs(path);
	if(!(ifs.is_open()))
	{
		std::cout<<"Can't read file "<<path<<std::endl;
		return;
	}
	std::string str;
	std::string index;
	std::stringstream ss;

	while(!ifs.eof())
	{
		str.clear();
		index.clear();
		ss.clear();

		std::getline(ifs,str);
		ss.str(str);
		ss>>index;

		if(!index.compare("v"))
		{
			double x,y,z;
			ss>>x>>y>>z;
			mVertices.push_back(Eigen::Vector3d(x,y,z));
		}
		else if(!index.compare("i"))
		{
			int i0,i1,i2,i3;
			ss>>i0>>i1>>i2>>i3;
			mTetrahedrons.push_back(Eigen::Vector4i(i0,i1,i2,i3));
		}
	}
	ifs.close();

	for(auto& v : mVertices)
		v = mT*v;
}