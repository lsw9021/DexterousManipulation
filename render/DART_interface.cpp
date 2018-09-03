#include "DART_interface.h"
using namespace dart::dynamics;
using namespace dart::simulation;

void
GUI::
DrawSkeleton(
	const dart::dynamics::SkeletonPtr& skel,
	const Eigen::Vector3d& color,
	bool box)
{
	for(int i=0;i<skel->getNumBodyNodes();i++)
	{
		auto bn = skel->getBodyNode(i);
		auto shapeNodes = bn->getShapeNodesWith<VisualAspect>();

		// for(int j=0;j<shapeNodes.size();j++){
		int j = (box? 1:0);
		if(shapeNodes.size() ==1)
			j=0;
		auto T = shapeNodes[j]->getTransform();
		DrawShape(T,shapeNodes[j]->getShape().get(),color);
		// }
	}
	// is_box = !is_box;
	if((skel->getPositions().norm()<1E-6))
	{
		for(int i =0;i<skel->getNumBodyNodes();i++)
		{
			Eigen::Isometry3d T;
			T = skel->getBodyNode(i)->getTransform();
			// std::cout<<T.linear()<<std::endl;
			// std::cout<<T.translation().transpose()<<std::endl;
			glPushMatrix();
			glMultMatrixd(T.data());
			glBegin(GL_LINES);
			glColor3f(1,0,0);
			glVertex3f(0,0,0);
			glVertex3f(0.1,0,0);

			glColor3f(0,1,0);
			glVertex3f(0,0,0);
			glVertex3f(0,0.1,0);

			glColor3f(0,0,1);
			glVertex3f(0,0,0);
			glVertex3f(0,0,0.1);

			glEnd();
			glPopMatrix();
		}
	}
	// for(int i =0;i<skel->getNumJoints();i++)
	// {
	// 	auto parent = skel->getJoint(i)->getParentBodyNode();
	// 	Eigen::Isometry3d T;
	// 	T.setIdentity();
	// 	if(parent!=nullptr)
	// 		T = parent->getTransform();
	// 	T = T*skel->getJoint(i)->getTransformFromParentBodyNode();
	// 	glPushMatrix();
	// 	glMultMatrixd(T.data());
	// 	glBegin(GL_LINES);
	// 	glColor3f(1,0,0);
	// 	glVertex3f(0,0,0);
	// 	glVertex3f(0.1,0,0);

	// 	glColor3f(0,1,0);
	// 	glVertex3f(0,0,0);
	// 	glVertex3f(0,0.1,0);

	// 	glColor3f(0,0,1);
	// 	glVertex3f(0,0,0);
	// 	glVertex3f(0,0,0.1);

	// 	glEnd();
	// 	glPopMatrix();
	// }
}


void
GUI::
DrawShape(const Eigen::Isometry3d& T,
	const dart::dynamics::Shape* shape,
	const Eigen::Vector3d& color)
{
	glEnable(GL_LIGHTING);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);
	glColor3f(color[0],color[1],color[2]);
	glPushMatrix();
	glMultMatrixd(T.data());
	if(shape->is<SphereShape>())
	{
		const auto* sphere = dynamic_cast<const SphereShape*>(shape);
		glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
		GUI::DrawSphere(sphere->getRadius());
		// glColor3f(0,0,0);
		// glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
		// GUI::DrawSphere(sphere->getRadius());
	}
	else if (shape->is<BoxShape>())
	{
		const auto* box = dynamic_cast<const BoxShape*>(shape);
		glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
    	GUI::DrawCube(box->getSize());
    	// GUI::DrawCube(Eigen::Vector3d(0.01,0.01,0.01));
	}
	else if(shape->is<MeshShape>())
	{
		auto* mesh = dynamic_cast<const MeshShape*>(shape);

		// for(int i =0;i<16;i++)
			// std::cout<<(*mesh->getMesh()->mRootNode->mTransformation)[i]<<" ";
    	GUI::DrawMesh(mesh->getScale(),mesh->getMesh());

	}

	glPopMatrix();

	// glDisable(GL_COLOR_MATERIAL);
}
void
GUI::
DrawMuscleWayPoints(const std::vector<AnchorPoint>& ap)
{
	std::vector<Eigen::Vector3d> point;

	for(int i=0;i<ap.size();i++)
		point.push_back(ap[i].first->getTransform()*ap[i].second);

	for(int i=0;i<ap.size()-1;i++)
		GUI::DrawLine(point[i],point[i+1],Eigen::Vector3d(0,0,0));
}