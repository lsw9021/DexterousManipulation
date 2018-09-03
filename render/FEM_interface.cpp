#include "FEM_interface.h"
#include "GL/glut.h"
using namespace FEM;
void
GUI::
DrawWorld(const std::shared_ptr<World>& world)
{
	const Eigen::VectorXd& X = world->GetPositions();
	const auto& cs = world->GetConstraints();

	for(const auto& c : cs)
		DrawConstraint(c, X);

}

void
GUI::
DrawConstraint(const std::shared_ptr<Cst>& c,const Eigen::VectorXd& x)
{
	glDisable(GL_LIGHTING);
	if(dynamic_cast<AttachmentCst*>(c.get()) != nullptr)
	{
		AttachmentCst* ac = dynamic_cast<AttachmentCst*>(c.get());	
		int i0 = ac->GetI0();
		const Eigen::Vector3d& p = ac->GetP();
		// glLineWidth(3.0);
		// GUI::DrawLine(x.block<3,1>(i0*3,0),p,Eigen::Vector3d(0.8,0,0));
		// glLineWidth(1.0);
	}
	else if(dynamic_cast<CorotateFEMCst*>(c.get()) != nullptr)
	{
		CorotateFEMCst* cc = dynamic_cast<CorotateFEMCst*>(c.get());
		int i0 = cc->GetI0();
		int i1 = cc->GetI1();
		int i2 = cc->GetI2();
		int i3 = cc->GetI3();

		const Eigen::Vector3d& p0 = x.block<3,1>(i0*3,0);
		const Eigen::Vector3d& p1 = x.block<3,1>(i1*3,0);
		const Eigen::Vector3d& p2 = x.block<3,1>(i2*3,0);
		const Eigen::Vector3d& p3 = x.block<3,1>(i3*3,0);

		GUI::DrawTetrahedron(p0,p1,p2,p3);
		glLineWidth(2.0);
		GUI::DrawLine(p0,p1,Eigen::Vector3d(0,0,0));
		GUI::DrawLine(p0,p2,Eigen::Vector3d(0,0,0));
		GUI::DrawLine(p0,p3,Eigen::Vector3d(0,0,0));
		GUI::DrawLine(p1,p2,Eigen::Vector3d(0,0,0));
		GUI::DrawLine(p1,p3,Eigen::Vector3d(0,0,0));
		GUI::DrawLine(p2,p3,Eigen::Vector3d(0,0,0));

		glLineWidth(1.0);
	}
	else if(dynamic_cast<LinearMuscleCst*>(c.get()) != nullptr)
	{
		LinearMuscleCst* cc = dynamic_cast<LinearMuscleCst*>(c.get());
		int i0 = cc->GetI0();
		int i1 = cc->GetI1();
		int i2 = cc->GetI2();
		int i3 = cc->GetI3();
		double a = cc->GetActivationLevel();
		const Eigen::Vector3d& p0 = x.block<3,1>(i0*3,0);
		const Eigen::Vector3d& p1 = x.block<3,1>(i1*3,0);
		const Eigen::Vector3d& p2 = x.block<3,1>(i2*3,0);
		const Eigen::Vector3d& p3 = x.block<3,1>(i3*3,0);

		GUI::DrawTetrahedron(p0,p1,p2,p3,Eigen::Vector3d(1.0,1.0-a,1.0-a));
		glLineWidth(2.0);
		GUI::DrawLine(p0,p1,Eigen::Vector3d(0,0,0));
		GUI::DrawLine(p0,p2,Eigen::Vector3d(0,0,0));
		GUI::DrawLine(p0,p3,Eigen::Vector3d(0,0,0));
		GUI::DrawLine(p1,p2,Eigen::Vector3d(0,0,0));
		GUI::DrawLine(p1,p3,Eigen::Vector3d(0,0,0));
		GUI::DrawLine(p2,p3,Eigen::Vector3d(0,0,0));

		glLineWidth(1.0);
	}
	else if(dynamic_cast<BarycentricAttachmentCst*>(c.get()) != nullptr)
	{
		BarycentricAttachmentCst* cc = dynamic_cast<BarycentricAttachmentCst*>(c.get());
		int i0 = cc->GetI0();
		int i1 = cc->GetI1();
		int i2 = cc->GetI2();
		int i3 = cc->GetI3();
		const Eigen::Vector3d& p0 = x.block<3,1>(i0*3,0);
		const Eigen::Vector3d& p1 = x.block<3,1>(i1*3,0);
		const Eigen::Vector3d& p2 = x.block<3,1>(i2*3,0);
		const Eigen::Vector3d& p3 = x.block<3,1>(i3*3,0);

		GUI::DrawTetrahedron(p0,p1,p2,p3,Eigen::Vector3d(1.0,0,0));
		glLineWidth(2.0);
		GUI::DrawLine(p0,p1,Eigen::Vector3d(0,0,0));
		GUI::DrawLine(p0,p2,Eigen::Vector3d(0,0,0));
		GUI::DrawLine(p0,p3,Eigen::Vector3d(0,0,0));
		GUI::DrawLine(p1,p2,Eigen::Vector3d(0,0,0));
		GUI::DrawLine(p1,p3,Eigen::Vector3d(0,0,0));
		GUI::DrawLine(p2,p3,Eigen::Vector3d(0,0,0));

		glLineWidth(1.0);
	}
	glEnable(GL_LIGHTING);
}