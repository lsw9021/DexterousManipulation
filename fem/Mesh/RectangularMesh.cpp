#include "RectangularMesh.h"

using namespace FEM;
RectangularMesh::
RectangularMesh(double _x,double _y, double _z, int _nx,int _ny, int _nz,const Eigen::Affine3d& _M)
	:mNx(_nx),mNy(_ny),mNz(_nz),mX(_x),mY(_y),mZ(_z),mT(_M)
{


	double dx = 0.5 * mX /(double)mNx;
	double dy = 0.5 * mY /(double)mNy;
	double dz = 0.5 * mZ /(double)mNz;

	double x,y,z;
	for(int i=0; i<2*mNx+1;i++)
	{	
		x = -0.5*mX + i*dx;
		for(int j=0; j<2*mNy+1;j++)
		{
			y = -0.5*mY + j*dy;
			for(int k=0; k<2*mNz+1;k++)
			{
				z = -0.5*mZ + k*dz;
				mVertices.push_back(Eigen::Vector3d(x,y,z));
			}
		}
	}

	for(auto& v : mVertices)
		v = mT*v;

	for(int i=0; i<mNx;i++)
	{	
		for(int j=0; j<mNy;j++)
		{
			for(int k=0; k<mNz;k++)
			{
				int index = (i*(2*mNy+1)*(2*mNz+1) + j*(2*mNz+1)+ k)*2;
				int cube[27] = {
					index,
					index+2*mNz+1,
					index+4*mNz+2,
					index+(2*mNz+1)*(2*mNy+1),
					index+(2*mNz+1)*(2*mNy+1)+2*mNz+1,
					index+(2*mNz+1)*(2*mNy+1)+4*mNz+2,
					index+(2*mNz+1)*(2*mNy+1)*2,
					index+(2*mNz+1)*(2*mNy+1)*2+2*mNz+1,
					index+(2*mNz+1)*(2*mNy+1)*2+4*mNz+2,
					index + 1,
					index+2*mNz+1 + 1,
					index+4*mNz+2 + 1,
					index+(2*mNz+1)*(2*mNy+1) + 1,
					index+(2*mNz+1)*(2*mNy+1)+2*mNz+1 + 1,
					index+(2*mNz+1)*(2*mNy+1)+4*mNz+2 + 1,
					index+(2*mNz+1)*(2*mNy+1)*2 + 1,
					index+(2*mNz+1)*(2*mNy+1)*2+2*mNz+1 + 1,
					index+(2*mNz+1)*(2*mNy+1)*2+4*mNz+2 + 1,
					index + 2,
					index+2*mNz+1 + 2,
					index+4*mNz+2 + 2,
					index+(2*mNz+1)*(2*mNy+1) + 2,
					index+(2*mNz+1)*(2*mNy+1)+2*mNz+1 + 2,
					index+(2*mNz+1)*(2*mNy+1)+4*mNz+2 + 2,
					index+(2*mNz+1)*(2*mNy+1)*2 + 2,
					index+(2*mNz+1)*(2*mNy+1)*2+2*mNz+1 + 2,
					index+(2*mNz+1)*(2*mNy+1)*2+4*mNz+2 + 2
				};

				int tets[][4] = {
					{0,1,3,9},
					{1,9,10,3},
					{3,9,12,10},
					{4,10,12,13},
					{1,3,4,10},
					{4,10,12,3},

					{1,4,5,10},
					{4,10,13,14},
					{4,14,10,5},
					{1,2,5,11},
					{10,11,14,5},
					{1,10,11,5},
				
					{3,6,7,15},
					{3,4,7,12},
					{7,12,15,3},
					{4,7,12,16},
					{4,12,13,16},
					{7,12,15,16},

					{4,13,16,14},
					{4,16,17,14},
					{4,7,16,17},

					{5,7,8,17},
					{5,7,17,4},
					{4,5,14,17},


					{12,21,19,22},
					{21,9,12,19},
					{9,18,19,21},
					{10,9,12,19},
					{12,13,10,22},
					{10,12,22,19},

					{10,13,14,22},
					{10,19,22,14},
					{11,14,19,10},
					{11,19,20,23},
					{23,22,19,14},
					{11,14,23,19},

					{15,21,24,25},
					{12,21,22,25},
					{12,15,21,25},
					{12,16,15,25},
					{12,13,16,22},
					{16,22,25,12},

					{14,16,17,25},
					{14,22,25,16},
					{13,14,16,22},
					{17,23,25,26},
					{14,17,22,23},
					{22,23,25,17}
				};
				for(int t = 0;t<48;t++)
				{
					int i0,i1,i2,i3;
					i0 = cube[tets[t][0]];
					i1 = cube[tets[t][1]];
					i2 = cube[tets[t][2]];
					i3 = cube[tets[t][3]];
					Eigen::Matrix3d Dm;
					Dm.block<3,1>(0,0) = mVertices[i1]-mVertices[i0];
					Dm.block<3,1>(0,1) = mVertices[i2]-mVertices[i0];
					Dm.block<3,1>(0,2) = mVertices[i3]-mVertices[i0];
					if(Dm.determinant()<0)
						mTetrahedrons.push_back(Eigen::Vector4i(i1,i0,i2,i3));
					else
						mTetrahedrons.push_back(Eigen::Vector4i(i0,i1,i2,i3));
				}
			}
		}
	}

}

std::shared_ptr<Mesh>
RectangularMesh::
Clone()
{
	return Create(mX,mY,mZ,mNx,mNy,mNz,mT);
}
std::shared_ptr<RectangularMesh>
RectangularMesh::
Create(double _x,double _y, double _z, int _nx,int _ny, int _nz,const Eigen::Affine3d& _M)
{
	auto rm = new RectangularMesh(_x,_y,_z,_nx,_ny,_nz,_M);

	return std::shared_ptr<RectangularMesh>(rm);
}
