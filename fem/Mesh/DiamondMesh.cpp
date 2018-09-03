#include "DiamondMesh.h"
#include <iostream>
#include <cassert>
using namespace FEM;

DiamondMesh::
DiamondMesh(double _x,double _y, double _z, int _nx,int _ny, int _nz,const Eigen::Affine3d& _M)
	:RectangularMesh(_x,_y,_z,_nx,_ny,_nz,Eigen::Affine3d::Identity())
{
	// mT = _M;
	// double dx = 0.5 * mX /(double)mNx;

	// mStartPointIndex = mVertices.size();
	// mEndPointIndex = mVertices.size() + 1;

	// double x_start = -0.5*mX - dx;
	// double x_end = -0.5*mX + (2*mNx+1)*dx;
	// mVertices.push_back(Eigen::Vector3d(x_start,0,0));
	// mVertices.push_back(Eigen::Vector3d(x_end,0,0));

	// int plus = 2*mNx*(2*mNy+1)*(2*mNz+1);
	// for(int i=0; i<mNy;i++)
	// {	
	// 	for(int j=0; j<mNz;j++)
	// 	{
			
	// 		int index = (i*(2*mNz+1)+ j)*2;
			
	// 		int tbt[9] = {
	// 			index,
	// 			index+1,
	// 			index+2,

	// 			index+(2*mNz+1),
	// 			index+(2*mNz+1)+1,
	// 			index+(2*mNz+1)+2,

	// 			index+(2*mNz+1)*2,
	// 			index+(2*mNz+1)*2+1,
	// 			index+(2*mNz+1)*2+2
	// 		};
	// 		mTetrahedrons.push_back(Eigen::Vector4i(mStartPointIndex,tbt[0],tbt[1],tbt[3]));
	// 		mTetrahedrons.push_back(Eigen::Vector4i(mStartPointIndex,tbt[1],tbt[3],tbt[4]));
	// 		mTetrahedrons.push_back(Eigen::Vector4i(mStartPointIndex,tbt[3],tbt[4],tbt[7]));
	// 		mTetrahedrons.push_back(Eigen::Vector4i(mStartPointIndex,tbt[3],tbt[6],tbt[7]));
	// 		mTetrahedrons.push_back(Eigen::Vector4i(mStartPointIndex,tbt[1],tbt[4],tbt[5]));
	// 		mTetrahedrons.push_back(Eigen::Vector4i(mStartPointIndex,tbt[1],tbt[2],tbt[5]));
	// 		mTetrahedrons.push_back(Eigen::Vector4i(mStartPointIndex,tbt[4],tbt[5],tbt[7]));
	// 		mTetrahedrons.push_back(Eigen::Vector4i(mStartPointIndex,tbt[7],tbt[5],tbt[8]));

	// 		mTetrahedrons.push_back(Eigen::Vector4i(mEndPointIndex,tbt[0]+plus,tbt[1]+plus,tbt[3]+plus));
	// 		mTetrahedrons.push_back(Eigen::Vector4i(mEndPointIndex,tbt[1]+plus,tbt[3]+plus,tbt[4]+plus));
	// 		mTetrahedrons.push_back(Eigen::Vector4i(mEndPointIndex,tbt[3]+plus,tbt[4]+plus,tbt[7]+plus));
	// 		mTetrahedrons.push_back(Eigen::Vector4i(mEndPointIndex,tbt[3]+plus,tbt[6]+plus,tbt[7]+plus));
	// 		mTetrahedrons.push_back(Eigen::Vector4i(mEndPointIndex,tbt[1]+plus,tbt[4]+plus,tbt[5]+plus));
	// 		mTetrahedrons.push_back(Eigen::Vector4i(mEndPointIndex,tbt[1]+plus,tbt[2]+plus,tbt[5]+plus));
	// 		mTetrahedrons.push_back(Eigen::Vector4i(mEndPointIndex,tbt[4]+plus,tbt[5]+plus,tbt[7]+plus));
	// 		mTetrahedrons.push_back(Eigen::Vector4i(mEndPointIndex,tbt[7]+plus,tbt[5]+plus,tbt[8]+plus));
	// 	}
	// }



	// for(int i =0;i<mTetrahedrons.size();i++)
	// {
	// 	auto& tets = mTetrahedrons[i];
	// 	int i0,i1,i2,i3;
	// 	i0 = tets[0];
	// 	i1 = tets[1];
	// 	i2 = tets[2];
	// 	i3 = tets[3];
	// 	Eigen::Matrix3d Dm;
	// 	Dm.block<3,1>(0,0) = mVertices[i1]-mVertices[i0];
	// 	Dm.block<3,1>(0,1) = mVertices[i2]-mVertices[i0];
	// 	Dm.block<3,1>(0,2) = mVertices[i3]-mVertices[i0];
	// 	if(Dm.determinant()<0)
	// 		mTetrahedrons[i] = Eigen::Vector4i(i1,i0,i2,i3);
		
	// }

	// double scale = 1.0/(x_end-x_start);
	// for(auto& v : mVertices)
	// 	v[0] = (scale)*v[0];

	
	mTetrahedrons.clear();
	mVertices.clear();

	mT = _M;
	mStartPointIndex = 9;
	mEndPointIndex = 8;
	mVertices = {
		Eigen::Vector3d(_x*0.25,-_y*0.5,_z*0.5),
		Eigen::Vector3d(_x*0.25,_y*0.5,_z*0.5),
		Eigen::Vector3d(_x*0.25,-_y*0.5,-_z*0.5),
		Eigen::Vector3d(_x*0.25,_y*0.5,-_z*0.5),
		Eigen::Vector3d(-_x*0.25,-_y*0.5,_z*0.5),
		Eigen::Vector3d(-_x*0.25,_y*0.5,_z*0.5),
		Eigen::Vector3d(-_x*0.25,-_y*0.5,-_z*0.5),
		Eigen::Vector3d(-_x*0.25,_y*0.5,-_z*0.5),
		Eigen::Vector3d(_x*0.5,0,0),
		Eigen::Vector3d(-_x*0.5,0,0)};

	mTetrahedrons = {
		Eigen::Vector4i(0,1,2,4),
		Eigen::Vector4i(1,2,3,7),
		Eigen::Vector4i(2,4,6,7),
		Eigen::Vector4i(1,4,5,7),
		Eigen::Vector4i(1,2,4,7),
		Eigen::Vector4i(0,1,2,8),
		Eigen::Vector4i(1,2,3,8),
		Eigen::Vector4i(4,5,7,9),
		Eigen::Vector4i(4,7,6,9)
	};


	for(auto& v : mVertices)
		v = mT*v;	
}

std::shared_ptr<Mesh>
DiamondMesh::
Clone()
{
	return Create(mX,mY,mZ,mNx,mNy,mNz,mT);
}
std::shared_ptr<DiamondMesh>
DiamondMesh::
Create(double _x,double _y, double _z, int _nx,int _ny, int _nz,const Eigen::Affine3d& _M)
{
	auto rm = new DiamondMesh(_x,_y,_z,_nx,_ny,_nz,_M);

	return std::shared_ptr<DiamondMesh>(rm);
}
