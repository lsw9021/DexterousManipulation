#include "IntegratedWorld.h"
#include "Record.h"
#include "MusculoSkeletalSystem.h"
#include "FSM.h"
#include "LLC.h"
#include <boost/filesystem.hpp>
#include <stdlib.h>
#include <memory>

int main(int argc,char** argv)
{
	std::string output_path("../output/");
	if(argc == 2)
		output_path = argv[1];
	if(output_path.back()!='/')
		output_path.append("/");
	if(boost::filesystem::exists(output_path))
		boost::filesystem::remove_all(output_path);
	boost::filesystem::create_directories(output_path);

	auto iw = std::make_shared<IntegratedWorld>();
	iw->Initialize();
	std::vector<int> V_list{
		3,3,3,3,3,
		3,3,3,4,4,
		4,4,4,4,4,
		4,4,4,5,5,
		5,5,5,5,5,
		5,5,5,5,5,
		5,5,5,5,5,
		5,5,5,5,5};
	std::shared_ptr<FSM> fsm = std::make_shared<FSM>(iw,V_list,5);
	std::shared_ptr<LowLevelController> llc = std::make_shared<LowLevelController>(iw);
	
	int dof = iw->GetMusculoSkeletalSystem()->GetSkeleton()->getNumDofs();
	Eigen::VectorXd p(dof),v(dof);
	p = iw->GetMusculoSkeletalSystem()->GetSkeleton()->getPositions();

	p.setZero();
	v.setZero();
	for(int i =0;i<2000;i++)
	{
		fsm->GetMotion(p,v);
		
		Eigen::VectorXd u = llc->Solve(p,v);
		iw->TimeStepping(u);

		iw->GetRecords().back()->Write(output_path+std::to_string(i));
	}
	
	system(("./render "+output_path).c_str());
	return 0;
}
