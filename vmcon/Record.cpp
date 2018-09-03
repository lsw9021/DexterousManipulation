#include "fem/fem.h"
#include "MusculoSkeletalSystem.h"
#include "Record.h"
#include <fstream>
#include <boost/filesystem.hpp>
using namespace dart::dynamics;
using namespace dart::simulation;
void
Record::
Set(const dart::simulation::WorldPtr& rigid_world,
	const std::shared_ptr<FEM::World>& soft_world,
	const std::shared_ptr<MusculoSkeletalSystem>& musculo_skeletal_system)
{
	t = rigid_world->getTime();
	for(int i =0;i<rigid_world->getNumSkeletons();i++)
		rigid_body_positions.insert(std::make_pair(rigid_world->getSkeleton(i)->getName(),rigid_world->getSkeleton(i)->getPositions()));

	for(int i =0;i<rigid_world->getNumSkeletons();i++)
		rigid_body_velocities.insert(std::make_pair(rigid_world->getSkeleton(i)->getName(),rigid_world->getSkeleton(i)->getVelocities()));

	soft_body_positions = soft_world->GetPositions();
	activation_levels = musculo_skeletal_system->GetActivationLevels();
}

void
Record::
Get(const dart::simulation::WorldPtr& rigid_world,
	const std::shared_ptr<FEM::World>& soft_world,
	const std::shared_ptr<MusculoSkeletalSystem>& musculo_skeletal_system)
{
	rigid_world->setTime(t);
	for(int i =0;i<rigid_world->getNumSkeletons();i++)
	{
		std::string key = rigid_world->getSkeleton(i)->getName();
		if(rigid_body_positions.find(key)!=rigid_body_positions.end())
			rigid_world->getSkeleton(i)->setPositions(rigid_body_positions.at(key));
		else
		{
			std::cout<<"No skeleton name "<<key<<std::endl;
			exit(0);
		}
	}
	for(int i =0;i<rigid_world->getNumSkeletons();i++)
	{
		std::string key = rigid_world->getSkeleton(i)->getName();
		if(rigid_body_velocities.find(key)!=rigid_body_velocities.end())
			rigid_world->getSkeleton(i)->setVelocities(rigid_body_velocities.at(key));
		else
		{
			std::cout<<"No skeleton name "<<key<<std::endl;
			exit(0);
		}
	}
	soft_world->SetTime(t);
	soft_world->SetPositions(soft_body_positions);
	int count = 0;
    for(auto& muscle : musculo_skeletal_system->GetMuscles())
    {
        muscle->SetActivationLevel(activation_levels[count]);
        muscle->origin->SetP(GetPoint(muscle->origin_way_points[0]));
        muscle->insertion->SetP(GetPoint(muscle->insertion_way_points[0]));
        count++;
    }
}

void
Record::
Write(const std::string& path)
{
	std::ofstream ofs(path);
	ofs<<"soft "<<soft_body_positions.transpose()<<std::endl;

	for(auto tup : rigid_body_positions)
	{
		ofs<<"rpos ";
		if(dart::math::isNan(tup.second))
		{
			ofs.close();
			boost::filesystem::remove(path);
			exit(0);
		}
		ofs<<tup.first<<" "<<tup.second.transpose()<<std::endl;	
	}
	for(auto tup : rigid_body_velocities)
	{
		ofs<<"rvel ";
		if(dart::math::isNan(tup.second))
		{
			ofs.close();
			boost::filesystem::remove(path);
			exit(0);
		}
		ofs<<tup.first<<" "<<tup.second.transpose()<<std::endl;	
	}

	ofs<<"time "<<t<<std::endl;
	ofs<<"act "<<activation_levels.transpose()<<std::endl;

	ofs.close();
}
void
Record::
Read(const std::string& path)
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

		Eigen::VectorXd eigen_vec;
		std::vector<double> vec;
		double val;
		if(!index.compare("soft"))
		{
			while(!ss.eof())
			{
				ss>>val;
				vec.push_back(val);
			}
			eigen_vec.resize(vec.size());
			for(int i=0;i<vec.size();i++)
				eigen_vec[i] = vec[i];
			soft_body_positions = eigen_vec;
		}
		else if(!index.compare("rpos"))
		{
			std::string name;
			ss>>name;
			while(!ss.eof())
			{
				ss>>val;
				vec.push_back(val);
			}
			eigen_vec.resize(vec.size());
			for(int i=0;i<vec.size();i++)
				eigen_vec[i] = vec[i];
			rigid_body_positions.insert(std::make_pair(name,eigen_vec));
		}
		else if(!index.compare("rvel"))
		{
			std::string name;
			ss>>name;
			while(!ss.eof())
			{
				ss>>val;
				vec.push_back(val);
			}
			eigen_vec.resize(vec.size());
			for(int i=0;i<vec.size();i++)
			{
				eigen_vec[i] = vec[i];
			}
			rigid_body_velocities.insert(std::make_pair(name,eigen_vec));

		}
		else if(!index.compare("act"))
		{
			while(!ss.eof())
			{
				ss>>val;
				vec.push_back(val);
			}
			eigen_vec.resize(vec.size());
			for(int i=0;i<vec.size();i++)
			{
				eigen_vec[i] = vec[i];
			}
			activation_levels = eigen_vec;
		}
		else if(!index.compare("time"))
		{
			ss>>val;
			t = val;
		}
		
	}
	ifs.close();
}