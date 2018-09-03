#include "MusculoSkeletalSystem.h"
#include "DART_helper.h"
#include <tinyxml.h>
#include <algorithm>
#include <fstream>
#include "Ball.h"
using namespace FEM;
using namespace dart::constraint;
using namespace dart::dynamics;
using namespace dart::simulation;
Eigen::Vector3d
GetPoint(const AnchorPoint& ap)
{
	return ap.first->getTransform()*ap.second;
}
int
Muscle::
GetNumForces()
{
	return origin_force_indices.size()+insertion_force_indices.size();
	// return 2;
}
Eigen::MatrixXd
Muscle::
GetJacobianTranspose()
{
	const auto& skel = origin_way_points[0].first->getSkeleton();
	Eigen::MatrixXd Jt(skel->getNumDofs(),3*GetNumForces());

	Jt.setZero();
	// std::cout<<Jt<<std::endl;
	// Jt.block(0,0,skel->getNumDofs(),3) = skel->getLinearJacobian(origin_way_points.back().first,origin_way_points.back().second).transpose();
	// Jt.block(0,3,skel->getNumDofs(),3) = skel->getLinearJacobian(insertion_way_points.back().first,insertion_way_points.back().second).transpose();
	int index = 0;
	for(int i : origin_force_indices){
		// std::cout<<i<<std::endl;
		// std::cout<< skel->getLinearJacobian(origin_way_points[i].first,origin_way_points[i].second).transpose()<<std::endl<<std::endl;
		Jt.block(0,index*3,skel->getNumDofs(),3) = skel->getLinearJacobian(origin_way_points[i].first,origin_way_points[i].second).transpose();
		// std::cout<<Jt<<std::endl<<std::endl;
		index++;
	}

	for(int i : insertion_force_indices){
		// std::cout<<i<<std::endl;
		// std::cout<< skel->getLinearJacobian(insertion_way_points[i].first,insertion_way_points[i].second).transpose()<<std::endl<<std::endl;
		Jt.block(0,index*3,skel->getNumDofs(),3) = skel->getLinearJacobian(insertion_way_points[i].first,insertion_way_points[i].second).transpose();
		// std::cout<<Jt<<std::endl<<std::endl;
		index++;
	}

	return Jt;
}
void
Muscle::
Initialize()
{
	int no = origin_way_points.size();
	int ni = insertion_way_points.size();

	if(!((origin_way_points[0].first == insertion_way_points[0].first) && (origin_way_points[0].first == origin_way_points[1].first)) )
		origin_force_indices.push_back(0);

	for(int i=1;i<no-1;i++)
	{
		if(!((origin_way_points[i].first == origin_way_points[i-1].first) && (origin_way_points[i].first == origin_way_points[i+1].first)) )
			origin_force_indices.push_back(i);
	}
	if(no>1)
	if( (origin_way_points[no-1].first != origin_way_points[no-2].first) )
	{
		origin_force_indices.push_back(no-1);
	}

	if(!((origin_way_points[0].first == insertion_way_points[0].first) && (insertion_way_points[0].first == insertion_way_points[1].first)) )
		insertion_force_indices.push_back(0);

	for(int i=1;i<ni-1;i++)
	{
		if(!((insertion_way_points[i].first == insertion_way_points[i-1].first) && (insertion_way_points[i].first == insertion_way_points[i+1].first)) )
			insertion_force_indices.push_back(i);
	}
	if(ni>1)
	if( (insertion_way_points[ni-1].first != insertion_way_points[ni-2].first) )
	{
		insertion_force_indices.push_back(ni-1);
	}
	// std::cout<<name<<std::endl;
	// for(int i : origin_force_indices)
	// 	std::cout<<i<<" ";
	// std::cout<<std::endl;
	// for(int i : insertion_force_indices)
	// 	std::cout<<i<<" ";
	// std::cout<<std::endl<<std::endl;
}

void
Muscle::
TransferForce(const Eigen::Vector3d& f_origin_tilda,const Eigen::Vector3d& f_insertion_tilda,Eigen::VectorXd& f)
{
	int no = origin_way_points.size();
	int ni = insertion_way_points.size();
	Eigen::VectorXd f_origin(no*3);
	Eigen::VectorXd f_insertion(ni*3);
	//origin
	f_origin.block<3,1>(0,0) = f_origin_tilda;
	Eigen::Vector3d u = (GetPoint(insertion_way_points[0])-GetPoint(origin_way_points[0])).normalized();
	for(int i=1;i<no;i++)
	{
		Eigen::Vector3d v = (GetPoint(origin_way_points[i-1])-GetPoint(origin_way_points[i])).normalized();
		double angle = acos(u.dot(v));
		Eigen::Vector3d axis = u.cross(v);
		axis.normalize();
		Eigen::AngleAxisd aa(angle,axis);
		f_origin.block<3,1>(i*3,0) = aa.toRotationMatrix()*f_origin_tilda;
	}
	//insertion
	f_insertion.block<3,1>(0,0) = f_insertion_tilda;
	u = -u;
	for(int i=1;i<ni;i++)
	{
		Eigen::Vector3d v = (GetPoint(insertion_way_points[i-1])-GetPoint(insertion_way_points[i])).normalized();
		double angle = acos(u.dot(v));
		Eigen::Vector3d axis = u.cross(v);
		axis.normalize();
		Eigen::AngleAxisd aa(angle,axis);
		f_insertion.block<3,1>(i*3,0) = aa.toRotationMatrix()*f_insertion_tilda;
	}

	f.resize(GetNumForces()*3);
	f.setZero();
	// f.block<3,1>(0,0) = f_origin.block<3,1>((no-1)*3,0);
	// f.block<3,1>(3,0) = f_insertion.block<3,1>((ni-1)*3,0);
	int index = 0;
	for(int i : origin_force_indices){
		auto before = std::find(origin_force_indices.begin(),origin_force_indices.end(),i-1);
		auto after = std::find(origin_force_indices.begin(),origin_force_indices.end(),i+1);
		// std::cout<<f.transpose()<<std::endl;
		if(i==0)
			f.block<3,1>(index*3,0) += f_origin.block<3,1>(i*3,0);
		// std::cout<<f.transpose()<<std::endl;
		
		if(before != origin_force_indices.end())
			f.block<3,1>(index*3,0) += f_origin.block<3,1>(i*3,0);
		// std::cout<<f.transpose()<<std::endl;

		if(after != origin_force_indices.end())
			f.block<3,1>(index*3,0) += -f_origin.block<3,1>((i+1)*3,0);
		// std::cout<<f.transpose()<<std::endl;

		index++;

	}

	for(int i : insertion_force_indices){
		auto before = std::find(insertion_force_indices.begin(),insertion_force_indices.end(),i-1);
		auto after = std::find(insertion_force_indices.begin(),insertion_force_indices.end(),i+1);
		// std::cout<<f.transpose()<<std::endl;

		if(i==0)
			f.block<3,1>(index*3,0) += f_insertion.block<3,1>(i*3,0);	
		// std::cout<<f.transpose()<<std::endl;
		
		if(before != insertion_force_indices.end())
			f.block<3,1>(index*3,0) += f_insertion.block<3,1>(i*3,0);
		// std::cout<<f.transpose()<<std::endl;
		if(after != insertion_force_indices.end())
			f.block<3,1>(index*3,0) += -f_insertion.block<3,1>((i+1)*3,0);
		// std::cout<<f.transpose()<<std::endl;
		index++;
	}
}

void
Muscle::
SetActivationLevel(double a)
{
	for(auto& lmc : muscle_csts)
		lmc->SetActivationLevel(a);

	activation_level = a;
}

MusculoSkeletalSystem::
MusculoSkeletalSystem()
	:mTendonStiffness(1E6),mMuscleStiffness(1E7),mYoungsModulus(5E6),mPoissonRatio(0.3)
{

}
void
MusculoSkeletalSystem::
AddMuscle(
	const std::string& name,
	const std::vector<AnchorPoint>& origin,
	const std::vector<AnchorPoint>& insertion,
	int origin_index,int insertion_index,
	const Eigen::Vector3d& fiber_direction,
	const MeshPtr& mesh)
{
	mMuscles.push_back(std::make_shared<Muscle>());
	auto& muscle = mMuscles.back();

	muscle->name = name;
	muscle->mesh = mesh;
	muscle->origin_way_points = origin;
	muscle->insertion_way_points = insertion;

	muscle->origin = AttachmentCst::Create(name+"_origin",mTendonStiffness,origin_index,GetPoint(origin[0]));
	muscle->insertion = AttachmentCst::Create(name+"_insertion",mTendonStiffness,insertion_index,GetPoint(insertion[0]));
	muscle->activation_level = 0.0;

	const auto& tetrahedrons = muscle->mesh->GetTetrahedrons();
	const auto& vertices = muscle->mesh->GetVertices();
	int tet_index = 0;
	for(const auto& tet: tetrahedrons)
	{
		int i0,i1,i2,i3;
		Eigen::Vector3d p0,p1,p2,p3;

		i0 = tet[0];
		i1 = tet[1];
		i2 = tet[2];
		i3 = tet[3];

		p0 = vertices[i0];
		p1 = vertices[i1];
		p2 = vertices[i2];
		p3 = vertices[i3];

		Eigen::Matrix3d Dm;

		Dm.block<3,1>(0,0) = p1 -p0;
		Dm.block<3,1>(0,1) = p2 -p0;
		Dm.block<3,1>(0,2) = p3 -p0;

		//Peventing inversion
		if(Dm.determinant()<0)
		{
			i2 = tet[3];
			i3 = tet[2];
			
			p2 = vertices[i2];
			p3 = vertices[i3];

			Dm.block<3,1>(0,1) = p2-p0;
			Dm.block<3,1>(0,2) = p3-p0;
		}

		muscle->muscle_csts.push_back(LinearMuscleCst::Create(name+"_muscle_"+std::to_string(tet_index),
			mMuscleStiffness,
			i0,i1,i2,i3,1.0/6.0*Dm.determinant(),Dm.inverse(),
			fiber_direction));
		muscle->csts.push_back(CorotateFEMCst::Create(name+"_element_"+std::to_string(tet_index),
			mYoungsModulus,
			mPoissonRatio,
			i0,i1,i2,i3,1.0/6.0*Dm.determinant(),Dm.inverse()));

		mAllMuscleConstraints.push_back(muscle->muscle_csts.back());
		tet_index++;
	}
	for(auto c: muscle->muscle_csts)
		muscle->csts.push_back(c);
	muscle->csts.push_back(muscle->origin);
	muscle->csts.push_back(muscle->insertion);
	muscle->Initialize();
}
void
MusculoSkeletalSystem::
Initialize(const std::shared_ptr<FEM::World>& soft_world,const dart::simulation::WorldPtr& rigid_world)
{
	for(int i =0;i<mMuscles.size();i++)
	{
		int offset = soft_world->GetNumVertices();
		auto& muscle = mMuscles[i];

		const auto& vertices = muscle->mesh->GetVertices();

		for(auto& c: muscle->csts)
			c->AddOffset(offset);
		Eigen::VectorXd v(vertices.size()*3);
		for(int i =0;i<vertices.size();i++)
			v.block<3,1>(i*3,0) = vertices[i];

		soft_world->AddBody(v,muscle->csts,1.0);
	}

	mActivationLevels.resize(mMuscles.size());
	mActivationLevels.setZero();

	rigid_world->addSkeleton(mSkeleton);
}
void
MusculoSkeletalSystem::
SetActivationLevels(const Eigen::VectorXd& a)
{
	mActivationLevels = a;
	for(int i =0;i<mMuscles.size();i++)
		mMuscles[i]->SetActivationLevel(a[i]);
}
void
MusculoSkeletalSystem::
TransformAttachmentPoints()
{
	for(auto& muscle : mMuscles)
	{
		auto& origin_way_points = muscle->origin_way_points;
		auto& insertion_way_points = muscle->insertion_way_points;

		Eigen::Vector3d po = GetPoint(origin_way_points[0]);
		Eigen::Vector3d pi = GetPoint(insertion_way_points[0]);

		muscle->origin->SetP(po);
		muscle->insertion->SetP(pi);
	}
}
void
MusculoSkeletalSystem::
ApplyForcesToSkeletons(const std::shared_ptr<FEM::World>& soft_world)
{
	// return;
	Eigen::VectorXd f = ComputeForce(soft_world);
	// std::cout<<f.transpose()<<std::endl;
	int index = 0;
	for(int i =0;i<mMuscles.size();i++)
	{
		auto& muscle = mMuscles[i];
		auto& origin_way_points = muscle->origin_way_points;
		auto& origin_force_indices = muscle->origin_force_indices;
		auto& insertion_way_points = muscle->insertion_way_points;
		auto& insertion_force_indices = muscle->insertion_force_indices;

		int n = muscle->GetNumForces();
		// origin_way_points.back().first->addExtForce(f.block<3,1>((index+0)*3,0),origin_way_points.back().second);
		// insertion_way_points.back().first->addExtForce(f.block<3,1>((index+1)*3,0),insertion_way_points.back().second);
		int no = origin_way_points.size();
		int count = 0;
		for(auto i : origin_force_indices)
		{
			// std::cout<<"origin "<<i<<"  "<<f.block<3,1>((index+count)*3,0).transpose()<<std::endl;
			origin_way_points[i].first->addExtForce(f.block<3,1>((index+count)*3,0),origin_way_points[i].second);
			count++;
		}
		for(auto i : insertion_force_indices)
		{
			// std::cout<<"insertion "<<i<<"  "<<f.block<3,1>((index+count)*3,0).transpose()<<std::endl;
			insertion_way_points[i].first->addExtForce(f.block<3,1>((index+count)*3,0),insertion_way_points[i].second);
			count++;
		}
		// std::cout<<std::endl;
		// for(int j=0;j<no;j++)
		// 	origin_way_points[j].first->addExtForce(f.block<3,1>((index+j)*3,0),origin_way_points[j].second);

		// for(int j=no;j<n;j++)
		// 	insertion_way_points[j-no].first->addExtForce(f.block<3,1>((index+j)*3,0),insertion_way_points[j-no].second);
		// origin_way_points.back().first->addExtForce(f.block<3,1>((index+0)*3,0),origin_way_points.back().second);
		// insertion_way_points.back().first->addExtForce(f.block<3,1>((index+1)*3,0),insertion_way_points.back().second);
		index +=n;
	}
}
void
MusculoSkeletalSystem::
ComputeForceDerivative(const FEM::WorldPtr& world,Eigen::SparseMatrix<double>& J)
{
	auto save_X = world->GetPositions();

	Eigen::VectorXd Ji(GetNumMuscleForces()*3);
	Eigen::VectorXd dg_da(save_X.rows()),dx_da(save_X.rows());
	Eigen::VectorXd temp_X(save_X.rows());
	Ji.setZero();
	dg_da.setZero();

#pragma omp parallel for
	for(int i=0;i<mAllMuscleConstraints.size();i++)
	{
		mAllMuscleConstraints[i]->Evaluatedgda(save_X);
	}
	for(int i=0;i<mAllMuscleConstraints.size();i++)
	{
		mAllMuscleConstraints[i]->Getdgda(dg_da);
	}

	world->Computedxda(dx_da,dg_da);

	temp_X = save_X + dx_da;
	world->SetPositions(temp_X);
	Ji = ComputeForce(world);
	world->SetPositions(save_X);
	Ji -= ComputeForce(world);




	std::vector<Eigen::Triplet<double>> j_triplets;

	j_triplets.reserve(GetNumMuscleForces()*3);
	int index = 0;
	for(int i =0;i<mMuscles.size();i++)
	{
		int n = mMuscles[i]->GetNumForces();
		for(int j=0;j<n;j++)
		{
			j_triplets.push_back(Eigen::Triplet<double>((index+j)*3+0,i,Ji[(index+j)*3+0]));
			j_triplets.push_back(Eigen::Triplet<double>((index+j)*3+1,i,Ji[(index+j)*3+1]));
			j_triplets.push_back(Eigen::Triplet<double>((index+j)*3+2,i,Ji[(index+j)*3+2]));
		}
		index +=n;
	}
	
	J.setFromTriplets(j_triplets.cbegin(), j_triplets.cend());
}
Eigen::VectorXd
MusculoSkeletalSystem::
ComputeForce(const FEM::WorldPtr& world)
{
	Eigen::VectorXd X = world->GetPositions();
	Eigen::VectorXd force_origin(X.rows()),force_insertion(X.rows());
	Eigen::VectorXd b(GetNumMuscleForces()*3);

	Eigen::Vector3d fo_tilda,fi_tilda;
	int index = 0;
	for(int i=0;i<mMuscles.size();i++)
	{
		auto& muscle = mMuscles[i];
		
		int n = muscle->GetNumForces();
		force_origin.setZero();
		force_insertion.setZero();

		muscle->origin->EvaluateGradient(X);
		muscle->insertion->EvaluateGradient(X);

		muscle->origin->GetGradient(force_origin);
		muscle->insertion->GetGradient(force_insertion);

		fo_tilda = force_origin.block<3,1>(muscle->origin->GetI0()*3,0);
		fi_tilda = force_insertion.block<3,1>(muscle->insertion->GetI0()*3,0);

		Eigen::VectorXd f;
		muscle->TransferForce(fo_tilda,fi_tilda,f);

		b.block(index*3,0,f.rows(),1) = f;
		index += n;
	}

	return b;
}
int
MusculoSkeletalSystem::
GetNumMuscleForces()
{
	int n=0;
	for(auto muscle : mMuscles)
		n += muscle->GetNumForces();
	return n;
}
void MakeMuscles(const std::string& path,std::shared_ptr<MusculoSkeletalSystem>& ms)
{
#ifdef USE_MUSCLE
	auto& skel = ms->GetSkeleton();
	TiXmlDocument doc;
    if(!doc.LoadFile(path))
    {
        std::cout<<"Cant open XML file : "<<path<<std::endl;
        return;
    }

    TiXmlElement* muscles = doc.FirstChildElement("Muscles");

    for(TiXmlElement* unit = muscles->FirstChildElement("unit");unit!=nullptr;unit = unit->NextSiblingElement("unit"))
    {
        TiXmlElement* ori = unit->FirstChildElement("origin");
        std::string name = (unit->Attribute("name"));
        std::vector<AnchorPoint> p_ori,p_ins;
       
        for(TiXmlElement* anc = ori->FirstChildElement("anchor");anc!=nullptr;anc = anc->NextSiblingElement("anchor"))   
        {
            std::string body_name = anc->Attribute("body");
            double x = std::stod(anc->Attribute("x"));
            double y = std::stod(anc->Attribute("y"));
            double z = std::stod(anc->Attribute("z"));

            auto T = skel->getBodyNode(body_name.c_str())->getShapeNodesWith<VisualAspect>()[0]->getRelativeTransform();
            auto T1 = skel->getBodyNode(body_name.c_str())->getTransform();
            
            Eigen::Vector3d body_coord(x,y,z);
            body_coord*=0.01;
            body_coord = T* body_coord;

            p_ori.push_back(AnchorPoint(skel->getBodyNode(body_name.c_str()),body_coord));
        }
        std::reverse(p_ori.begin(),p_ori.end());
        TiXmlElement* ins = unit->FirstChildElement("insertion");
        for(TiXmlElement* anc = ins->FirstChildElement("anchor");anc!=nullptr;anc = anc->NextSiblingElement("anchor"))   
        {
            std::string body_name = anc->Attribute("body");
            double x = std::stod(anc->Attribute("x"));
            double y = std::stod(anc->Attribute("y"));
            double z = std::stod(anc->Attribute("z"));

            auto T = skel->getBodyNode(body_name.c_str())->getShapeNodesWith<VisualAspect>()[0]->getRelativeTransform();
            auto T1 = skel->getBodyNode(body_name.c_str())->getTransform();
            
            Eigen::Vector3d body_coord(x,y,z);
            body_coord*=0.01;
            body_coord = T* body_coord;

            p_ins.push_back(AnchorPoint(skel->getBodyNode(body_name.c_str()),body_coord));
        }

        Eigen::Vector3d muscle_start,muscle_end;

        muscle_start = GetPoint(p_ori[0]);
        muscle_end = GetPoint(p_ins[0]);

        double len = (muscle_start - muscle_end).norm();
        Eigen::Vector3d unit_dir = (muscle_start - muscle_end).normalized();

        Eigen::Vector3d axis = Eigen::Vector3d::UnitX().cross(unit_dir);
        double cos_angle = unit_dir[0];
        double sin_angle = axis.norm();

        double angle = atan2(sin_angle,cos_angle);
        TiXmlElement* mesh_element = unit->FirstChildElement("mesh");

        Eigen::Affine3d T = Eigen::Affine3d::Identity();
        T.translation() = 0.5*(muscle_start + muscle_end);
        T.linear() = len*(Eigen::AngleAxisd(angle,axis.normalized()).matrix());
        int nx = std::stoi(mesh_element->Attribute("nx"));
        int ny = std::stoi(mesh_element->Attribute("ny"));
        double ratio = std::stod(mesh_element->Attribute("ratio"));
        auto dm = DiamondMesh::Create(1.0,(double)ny/(double)nx*ratio,(double)ny/(double)nx*ratio,nx,ny,ny,T);
        int i_ori = dm->GetEndingPointIndex();
        int i_ins = dm->GetStartingPointIndex();

        ms->AddMuscle(name,p_ori,p_ins,i_ori,i_ins,unit_dir,dm);
        
        
    }
#endif
}
void MakeSkeleton(std::shared_ptr<MusculoSkeletalSystem>& ms)
{
	ms->GetSkeleton() = Skeleton::create("HUMAN");
	auto& skel = ms->GetSkeleton();

	MakeRootBody(skel,"Torso",
		Eigen::Vector3d(0.03,0.6,0.03),
		Eigen::Vector3d(0,-0.3,0),
		JOINT_TYPE::BALL_AND_SOCKET,10);


	MakeBody(skel,skel->getBodyNode("Torso"),"NeckR",
		Eigen::Vector3d(0.2,0.03,0.03),
		Eigen::Vector3d(0.0,0.15,0),
		Eigen::Vector3d(0.1,0.0,0.0),JOINT_TYPE::REVOLUTE,5);

	MakeBody(skel,skel->getBodyNode("Torso"),"NeckL",
		Eigen::Vector3d(0.2,0.03,0.03),
		Eigen::Vector3d(0.0,0.15,0),
		Eigen::Vector3d(-0.1,0.0,0),JOINT_TYPE::REVOLUTE,5);

	MakeBody(skel,skel->getBodyNode("NeckR"),"ShoulderR",
		Eigen::Vector3d(0.3,0.03,0.03),
		Eigen::Vector3d(-0.08,-0.03,-0.03),
		Eigen::Vector3d(0.15,0.0,0),JOINT_TYPE::BALL_AND_SOCKET,5);

	MakeBody(skel,skel->getBodyNode("NeckL"),"ShoulderL",
		Eigen::Vector3d(0.3,0.03,0.03),
		Eigen::Vector3d(0.08,-0.03,-0.03),
		Eigen::Vector3d(-0.15,0.0,0),JOINT_TYPE::BALL_AND_SOCKET,5);

	MakeBody(skel,skel->getBodyNode("ShoulderR"),"ElbowR",
		Eigen::Vector3d(0.3,0.03,0.03),
		Eigen::Vector3d(-0.185,0.0,0.02),
		Eigen::Vector3d(0.145,0.0,0),JOINT_TYPE::REVOLUTE,5);
	
	MakeBody(skel,skel->getBodyNode("ShoulderL"),"ElbowL",
		Eigen::Vector3d(0.3,0.03,0.03),
		Eigen::Vector3d(0.185,0.0,0.02),
		Eigen::Vector3d(-0.145,0.0,0),JOINT_TYPE::REVOLUTE,5);

	MakeBody(skel,skel->getBodyNode("ElbowR"),"HandR",
		Eigen::Vector3d(0.07,0.07,0.07),
		Eigen::Vector3d(-0.13,0.05,0.0),
		Eigen::Vector3d(0,0,0),
		JOINT_TYPE::BALL_AND_SOCKET,
		3);

	MakeBody(skel,skel->getBodyNode("ElbowL"),"HandL",
		Eigen::Vector3d(0.07,0.07,0.07),
		Eigen::Vector3d(0.13,0.05,0.0),
		Eigen::Vector3d(0,0,0),
		JOINT_TYPE::BALL_AND_SOCKET,
		3);
	MakeBody(skel,skel->getBodyNode("Torso"),"Head",
		Eigen::Vector3d(0.07,0.07,0.07),
		Eigen::Vector3d(0,0.33,0),
		Eigen::Vector3d(0,0,0),
		JOINT_TYPE::WELD,
		10);
	Eigen::VectorXd pos = skel->getPositions();
	int dt = skel->getJoint(0)->getNumDofs();
	int dn = skel->getJoint(1)->getNumDofs();

	//Torso Joint
	pos[0] = 0.0;
	pos[1] = 0.0;
	pos[2] = 0.0;
	pos[3] = 0.0;
	pos[4] = 0.0;
	pos[5] = 0.0;
	//NeckR Joint
	pos[dt+0] = 0.0;
	pos[dt+1] = 0.0;
	pos[dt+2] = 0.0;
	//NeckL Joint
	pos[dt+dn+0] = -0.0;
	pos[dt+dn+1] = -0.0;
	pos[dt+dn+2] = -0.0;
	//ShoulderR Joint
	pos[dt+dn*2+0] =0.3;
	pos[dt+dn*2+1] =0.0;
	pos[dt+dn*2+2] =1.2;
	//ShoulderL Joint
	pos[dt+dn*2+3] =0.3;
	pos[dt+dn*2+4] =-0.0;
	pos[dt+dn*2+5] =-1.2;
	//ElbowR Joint
	pos[dt+dn*2+6] = 1.7;
	//ElbowL Joint
	pos[dt+dn*2+7] = -1.7;
	//HandR Joint
	pos[dt+dn*2+8] = 0.0;
	pos[dt+dn*2+9] = 0.0;
	pos[dt+dn*2+10] = 0.0;
	//HandL Joint
	pos[dt+dn*2+11] = 0.0;
	pos[dt+dn*2+12] = 0.0;
	pos[dt+dn*2+13] = 0.0;

	//Torso joint
	skel->getDof(0)->setPositionLimits(-0.0,0.0);
	skel->getDof(1)->setPositionLimits(-0.0,0.0);
	skel->getDof(2)->setPositionLimits(-0.0,0.0);
	//NeckR Joint
	
	skel->getDof(dt+0)->setPositionLimits(-0.0,0.0);
	//NeckL Joint
	skel->getDof(dt+dn+0)->setPositionLimits(-0.0,0.0);
	//ShoulderR Joint
	//ShoulderL Joint
	//ElbowR Joint
	skel->getDof(dt+dn*2+6)->setPositionLimits(0.1,2.2);
	//ElbowL Joint
	skel->getDof(dt+dn*2+7)->setPositionLimits(-2.2,0.1); 
	//HandR Joint
	skel->getDof(dt+dn*2+8)->setPositionLimits(-1.0,1.0);
	skel->getDof(dt+dn*2+9)->setPositionLimits(-1.0,1.0);
	skel->getDof(dt+dn*2+10)->setPositionLimits(-0.8,0.8);
	//HandL Joint
	skel->getDof(dt+dn*2+11)->setPositionLimits(-1.0,1.0);
	skel->getDof(dt+dn*2+12)->setPositionLimits(-1.0,1.0);
	skel->getDof(dt+dn*2+13)->setPositionLimits(-0.8,0.8);

	for(int i =0;i<skel->getNumDofs();i++){
		skel->getDof(i)->getJoint()->setPositionLimitEnforced(true);
	}
	// skel->setMobile(false);
	for(int i=0;i<skel->getNumBodyNodes();i++)
		skel->getBodyNode(i)->setCollidable(false);

	skel->setPositions(pos);
	skel->computeForwardKinematics(true,false,false);

}

void MakeBalls(dart::simulation::WorldPtr& world,const std::shared_ptr<MusculoSkeletalSystem>& ms,std::vector<std::shared_ptr<Ball>>& ball,int num)
{
	for(int i =0;i<num;i++)
	{
		SkeletonPtr skel = Skeleton::create("ball_"+std::to_string(i));

		bool is_left_hand = i%2;
		if(is_left_hand)
		{
			auto* abn =ms->GetSkeleton()->getBodyNode("HandL");
			Eigen::Vector3d loc = abn->getTransform().translation();
			loc += Eigen::Vector3d(-0.05,0.015,0.10-0.009*i);
			MakeBall(skel,loc,0.036,0.13);

			for(int i =0;i<skel->getNumDofs();i++){
				//skel->getDof(i)->setDampingCoefficient(0.1);
			}
			ball.push_back(std::make_shared<Ball>(skel));
			world->addSkeleton(skel);
			ball.back()->Attach(world,abn);
		}
		else
		{
			auto* abn =ms->GetSkeleton()->getBodyNode("HandR");
			Eigen::Vector3d loc = abn->getTransform().translation();
			loc += Eigen::Vector3d(0.05,0.015,0.10-0.009*i);
			MakeBall(skel,loc,0.036,0.13);
			for(int i =0;i<skel->getNumDofs();i++){
				//skel->getDof(i)->setDampingCoefficient(0.1);
			}
			ball.push_back(std::make_shared<Ball>(skel));
			world->addSkeleton(skel);
			ball.back()->Attach(world,abn);	
		}
	}

}
