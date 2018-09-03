#include "World.h"
#include <chrono>
using namespace FEM;
#define BIG_EPSILON 1E-4
#define EPSILON 1E-6
World::
World(	IntegrationMethod im,
		double time_step,
		int max_iteration,
		const Eigen::Vector3d& gravity,
		double damping_coeff)
	:mIntegrationMethod(im),
	mTimeStep(time_step),
	mMaxIteration(max_iteration),
	mGravity(gravity),
	mDampingCoefficient(damping_coeff),
	mNumVertices(0),
	mConstraintDofs(0),
	mIsInitialized(false)
{

}
std::shared_ptr<World>
World::
Clone()
{
	auto new_world = Create();

	new_world->mIsInitialized = mIsInitialized;
	new_world->mNumVertices = mNumVertices;
	new_world->mConstraintDofs = mConstraintDofs;
	new_world->mMaxIteration = mMaxIteration;

	new_world->mTimeStep = mTimeStep;
	new_world->mTime = mTime;
	new_world->mDampingCoefficient = mDampingCoefficient;
	new_world->mGravity = mGravity;

	new_world->mIntegrationMethod = mIntegrationMethod;

	new_world->mUnitMass = mUnitMass;
	new_world->mConstraints = mConstraints;

	for(int i=0;i<mConstraints.size();i++)
		new_world->mConstraints[i] = mConstraints[i]->Clone();

	new_world->mPositions = mPositions;
	new_world->mVelocities = mVelocities;
	new_world->mExternalForces = mExternalForces;

	new_world->mMassMatrix = mMassMatrix;
	new_world->mInvMassMatrix = mInvMassMatrix;
	new_world->mIdentityMatrix = mIdentityMatrix;
	
	new_world->mq = mq;
	new_world->mJ = mJ;
	new_world->mL = mL;

	if(mIntegrationMethod == PROJECTIVE_DYNAMICS)
	{
		Eigen::SparseMatrix<double> H2ML = (1.0/(mTimeStep*mTimeStep))*mMassMatrix+mL;
		FactorizeLDLT(H2ML,new_world->mSolver);
	}
	else if(mIntegrationMethod == PROJECTIVE_QUASI_STATIC)
		FactorizeLDLT(mL,new_world->mSolver);

	return new_world;
}
std::shared_ptr<World>
World::
Create(	IntegrationMethod im,
		double time_step,
		int max_iteration,
		const Eigen::Vector3d& gravity,
		double damping_coeff)
{
	auto w = new World(im,time_step,max_iteration,gravity,damping_coeff);
	return std::shared_ptr<World>(w);
}
void
World::
Initialize()
{
	mTime = 0.0;
	
	mVelocities.resize(mNumVertices*3);	
	mVelocities.setZero();

	mExternalForces.resize(mNumVertices*3);
	mExternalForces.setZero();

	mMassMatrix.resize(mNumVertices*3,mNumVertices*3);
	mInvMassMatrix.resize(mNumVertices*3,mNumVertices*3);
	mIdentityMatrix.resize(mNumVertices*3,mNumVertices*3);

	std::vector<Eigen::Triplet<double>> i_triplets;
	std::vector<Eigen::Triplet<double>> m_triplets;
	std::vector<Eigen::Triplet<double>> inv_m_triplets;

	for(int i =0;i<mNumVertices;i++)
	{
		m_triplets.push_back(Eigen::Triplet<double>(i*3+0,i*3+0,mUnitMass[i]));
		m_triplets.push_back(Eigen::Triplet<double>(i*3+1,i*3+1,mUnitMass[i]));
		m_triplets.push_back(Eigen::Triplet<double>(i*3+2,i*3+2,mUnitMass[i]));

		inv_m_triplets.push_back(Eigen::Triplet<double>(i*3+0,i*3+0,1.0/mUnitMass[i]));
		inv_m_triplets.push_back(Eigen::Triplet<double>(i*3+1,i*3+1,1.0/mUnitMass[i]));
		inv_m_triplets.push_back(Eigen::Triplet<double>(i*3+2,i*3+2,1.0/mUnitMass[i]));

		i_triplets.push_back(Eigen::Triplet<double>(i*3+0,i*3+0,1.0));
		i_triplets.push_back(Eigen::Triplet<double>(i*3+1,i*3+1,1.0));
		i_triplets.push_back(Eigen::Triplet<double>(i*3+2,i*3+2,1.0));
	}

	mMassMatrix.setFromTriplets(m_triplets.cbegin(), m_triplets.cend());
	mInvMassMatrix.setFromTriplets(inv_m_triplets.cbegin(), inv_m_triplets.cend());
	mIdentityMatrix.setFromTriplets(i_triplets.cbegin(), i_triplets.cend());

	mq.resize(mNumVertices*3);
	mq.setZero();

	if( mIntegrationMethod == IntegrationMethod::PROJECTIVE_DYNAMICS ||
		mIntegrationMethod == IntegrationMethod::PROJECTIVE_QUASI_STATIC)
		Precompute();

	mIsInitialized = true;

	std::cout<<"Total degree of freedom : "<<mPositions.rows()<<std::endl;
	std::cout<<"Total constraints : "<<mConstraints.size()<<std::endl;
}
void
World::
TimeStepping(bool integrate)
{
	if(!mIsInitialized)
	{
		std::cout<<"Engine is not initialized."<<std::endl;
		return;
	}
	Eigen::VectorXd x_next(mNumVertices*3);
	mExternalForces.setZero();
	for(int i=0;i<mNumVertices;i++)
		mExternalForces.block<3,1>(i*3,0) = mGravity;

	mq = mPositions + mTimeStep*mVelocities + (mTimeStep*mTimeStep)*(mInvMassMatrix*mExternalForces);
	switch(mIntegrationMethod)
	{
	case NEWTON_METHOD:
		x_next = IntegrateNewtonMethod();
		break;
	case QUASI_STATIC:
		x_next = IntegrateQuasiStatic();
		break;
	case PROJECTIVE_DYNAMICS:
	
		x_next = IntegrateProjectiveDynamics();
		break;
	case PROJECTIVE_QUASI_STATIC:
		x_next = IntegrateProjectiveQuasiStatic();
		break;
	default:
	return;
	};
	
	if(integrate)
	{
		IntegratePositionsAndVelocities(x_next);
		mTime += mTimeStep;
	}
	else
		mPositions = x_next;
}
void
World::
AddBody(const Eigen::VectorXd& x0, const std::vector<std::shared_ptr<Cst>>& constraints,double m)
{
	if(mIsInitialized)
	{
		std::cout<<"Add Body before initializing engine."<<std::endl;
		return;
	}
	int nv = x0.rows()/3;
	mNumVertices += nv;

	Eigen::VectorXd prev_x(mPositions);
	mPositions.resize(mNumVertices*3);
	mPositions.head(prev_x.rows()) = prev_x;
	mPositions.tail(x0.rows()) = x0;

	mConstraints.insert(mConstraints.end(),constraints.begin(),constraints.end());
	double unit_mass = m/((double)nv);
	for(int i=0;i<nv;i++)
		mUnitMass.push_back(unit_mass);
}
void
World::
AddConstraint(const std::shared_ptr<Cst>& c)
{
	mConstraints.push_back(c);
	if((mIntegrationMethod == IntegrationMethod::PROJECTIVE_DYNAMICS||
		mIntegrationMethod == IntegrationMethod::PROJECTIVE_QUASI_STATIC)&&
		mIsInitialized)
	{
		mConstraintDofs = 0;
		Precompute();
	}
}
void
World::
RemoveConstraint(const std::shared_ptr<Cst>& c)
{
	bool isRemoved = false;
	for(int i =0;i<mConstraints.size();i++)
	{
		if(mConstraints[i]==c)
		{
			mConstraints.erase(mConstraints.begin()+i);
			isRemoved = true;
			break;
		}
	}

	if((mIntegrationMethod == IntegrationMethod::PROJECTIVE_DYNAMICS||
		mIntegrationMethod == IntegrationMethod::PROJECTIVE_QUASI_STATIC)&&
		mIsInitialized)
	{
		Precompute();
	}
}
void
World::
Computedxda(Eigen::VectorXd& dx_da,const Eigen::VectorXd& dg_da)
{
	if(mIntegrationMethod==PROJECTIVE_QUASI_STATIC)
		dx_da = -mSolver.solve(dg_da);
	else
	{
		Eigen::SparseMatrix<double> H(mNumVertices*3,mNumVertices*3);
		EvaluateConstraintsHessian(mPositions,H);
		FactorizeLDLT(H,mSolver);
		dx_da = -mSolver.solve(dg_da);
	}
}
Eigen::VectorXd
World::
IntegrateNewtonMethod()
{
	Eigen::VectorXd x_next(mNumVertices*3);

	x_next = mq;

	Eigen::VectorXd g_k(mNumVertices*3);
	Eigen::SparseMatrix<double> H_k(mNumVertices*3,mNumVertices*3);
	for(int k=0;k<mMaxIteration;k++)
	{
		EvaluateGradient(x_next,g_k);

		if(g_k.squaredNorm()<BIG_EPSILON)
			break;

		EvaluateHessian(x_next,H_k);
		FactorizeLDLT(H_k,mSolver);
		Eigen::VectorXd d = -mSolver.solve(g_k);

		double alpha = ComputeStepSize(x_next,g_k,d);

		x_next += alpha*d;
	}
	InversionFree(x_next);
	return x_next;
}
Eigen::VectorXd
World::
IntegrateQuasiStatic()
{
	Eigen::VectorXd x_next(mNumVertices*3);


	x_next = mPositions;
	Eigen::VectorXd g_k(mNumVertices*3);
	Eigen::SparseMatrix<double> H_k(mNumVertices*3,mNumVertices*3);
	for(int k=0;k<mMaxIteration;k++)
	{
		EvaluateConstraintsGradient(x_next,g_k);

		if(g_k.squaredNorm()<BIG_EPSILON)
			break;

		EvaluateConstraintsHessian(x_next,H_k);
		FactorizeLDLT(H_k,mSolver);
		Eigen::VectorXd d = -mSolver.solve(g_k);

		double alpha = ComputeStepSize(x_next,g_k,d);

		x_next += alpha*d;
	}
	InversionFree(x_next);
	return x_next;
}
Eigen::VectorXd
World::
IntegrateProjectiveDynamics()
{
	Eigen::VectorXd x_next(mNumVertices*3);
	Eigen::VectorXd b(mNumVertices*3);
	Eigen::VectorXd d(mConstraintDofs*3);

	b = (1.0/(mTimeStep*mTimeStep))*mMassMatrix*mq;

	x_next = mq;

	for(int k=0;k<mMaxIteration;k++)
	{

		EvaluateDVector(x_next,d);
		
		x_next = mSolver.solve(b+mJ*d);
		
	}
	InversionFree(x_next);
	return x_next;
}
Eigen::VectorXd
World::
IntegrateProjectiveQuasiStatic()
{
	Eigen::VectorXd x_next(mNumVertices*3);
	Eigen::VectorXd d(mConstraintDofs*3);

	x_next = mPositions;

	for(int k=0;k<mMaxIteration;k++)
	{
		EvaluateDVector(x_next,d);
		x_next = mSolver.solve(mJ*d);
	}
	InversionFree(x_next);
	return x_next;
}
void
World::
FactorizeLLT(const Eigen::SparseMatrix<double>& A, Eigen::SimplicialLLT<Eigen::SparseMatrix<double>>& lltSolver)
{
	Eigen::SparseMatrix<double> A_prime = A;
	lltSolver.analyzePattern(A_prime);
	lltSolver.factorize(A_prime);
	double damping = 1E-6;
	bool success = true;
	while (lltSolver.info() != Eigen::Success)
	{
		damping *= 10;
		A_prime = A + damping*mIdentityMatrix;
		lltSolver.factorize(A_prime);
		success = false;
	}
	if (!success)
		std::cout << "factorize failure (damping : " << damping<<" )"<<std::endl;
}
void
World::
FactorizeLDLT(const Eigen::SparseMatrix<double>& A, Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>& ldltSolver)
{
	Eigen::SparseMatrix<double> A_prime = A;
	ldltSolver.analyzePattern(A_prime);
	ldltSolver.factorize(A_prime);
	double damping = 1E-6;
	bool success = true;
	while (ldltSolver.info() != Eigen::Success)
	{
		damping *= 10;
		A_prime = A + damping*mIdentityMatrix;
		ldltSolver.factorize(A_prime);
		success = false;
	}
	if (!success)
		std::cout << "factorize failure (damping : " << damping<<" )"<<std::endl;
}
double
World::
EvaluateEnergy(const Eigen::VectorXd& x)
{
	Eigen::VectorXd x_q = x - mq;

	double val = EvaluateConstraintsEnergy(x);

	if(mIntegrationMethod == NEWTON_METHOD)
		val += 0.5*(1.0/(mTimeStep*mTimeStep))*(x_q.dot(mMassMatrix*x_q));

	return val;
}
void
World::
EvaluateGradient(const Eigen::VectorXd& x,Eigen::VectorXd& g)
{
	EvaluateConstraintsGradient(x,g);
	
	if(mIntegrationMethod == NEWTON_METHOD)
		g += (1.0/(mTimeStep*mTimeStep))*mMassMatrix*(x - mq);
}
void
World::
EvaluateHessian(const Eigen::VectorXd& x,Eigen::SparseMatrix<double>& H)
{
	EvaluateConstraintsHessian(x,H);

	if(mIntegrationMethod == NEWTON_METHOD)
		H = H + (1.0/(mTimeStep*mTimeStep))*mMassMatrix;
}

double
World::
EvaluateConstraintsEnergy(const Eigen::VectorXd& x)
{
	double energy=0;
	
#pragma omp parallel for
	for(int i=0;i<mConstraints.size();i++)
	{
		mConstraints[i]->EvaluatePotentialEnergy(x);
	}

	for(auto& c : mConstraints)
	{
		c->GetPotentialEnergy(energy);
	}


	return energy;
}
void
World::
EvaluateConstraintsGradient(const Eigen::VectorXd& x,Eigen::VectorXd& g)
{
	g.resize(mNumVertices*3);
	g.setZero();

#pragma omp parallel for
	for(int i=0;i<mConstraints.size();i++)
	{
		mConstraints[i]->EvaluateGradient(x);
	}

	for(auto& c : mConstraints)
	{
		c->GetGradient(g);
	}

}
void
World::
EvaluateConstraintsHessian(const Eigen::VectorXd& x,Eigen::SparseMatrix<double>& H)
{
	H.resize(mNumVertices*3,mNumVertices*3);
	std::vector<Eigen::Triplet<double>> h_triplets;

#pragma omp parallel for
	for(int i=0;i<mConstraints.size();i++)
	{
		mConstraints[i]->EvaluateHessian(x);
	}

	for(auto& c : mConstraints)
	{
		c->GetHessian(h_triplets);
	}

	H.setFromTriplets(h_triplets.cbegin(),h_triplets.cend());
}
double
World::
ComputeStepSize(const Eigen::VectorXd& x, const Eigen::VectorXd& g,const Eigen::VectorXd& d)
{
	double alpha = 1.0;
	double c1 = 0.00;
	double c2 = 0.1;
	double f_x,f_x_next;
	Eigen::VectorXd x_next;

	f_x = EvaluateEnergy(x);

	double g_dot_d = g.dot(d);

	for(int i=0;i<16;i++)
	{
		x_next = x + alpha*d;

		f_x_next = EvaluateEnergy(x_next);

		if(f_x + alpha*c1*g_dot_d>=f_x_next)
			break;

		alpha *= c2;
	}

	return alpha;
}
void
World::
EvaluateDVector(const Eigen::VectorXd& x,Eigen::VectorXd& d)
{
	d.resize(mConstraintDofs*3);

	int n = mConstraints.size();
#pragma omp parallel for
	for(int i=0;i<n;i++)
	{
		mConstraints[i]->EvaluateDVector(x);
	}
	int index = 0;
	for(auto& c : mConstraints)
	{
		c->GetDVector(index,d);
	}
}
void
World::
EvaluateJMatrix(Eigen::SparseMatrix<double>& J)
{
	J.resize(mNumVertices*3,mConstraintDofs*3);
	std::vector<Eigen::Triplet<double>> J_triplets;

	int index = 0;
	for(auto& c : mConstraints)
	{
		c->EvaluateJMatrix(index,J_triplets);
	}

	J.setFromTriplets(J_triplets.cbegin(),J_triplets.cend());
}
void
World::
EvaluateLMatrix(Eigen::SparseMatrix<double>& L)
{
	L.resize(mNumVertices*3,mNumVertices*3);
	std::vector<Eigen::Triplet<double>> L_triplets;

	for(auto& c: mConstraints)
	{
		c->EvaluateLMatrix(L_triplets);
	}

	L.setFromTriplets(L_triplets.cbegin(),L_triplets.cend());
}
void
World::
Precompute()
{
	mConstraintDofs = 0;

	//For computing constraint's dofs
	Eigen::VectorXd d_temp(3*mConstraints.size()*6);
	int index = 0;
	for(auto& c : mConstraints)
	{
		c->GetDVector(index,d_temp);
	}
	mConstraintDofs = index;

	EvaluateLMatrix(mL);
	EvaluateJMatrix(mJ);

	Eigen::SparseMatrix<double> H2ML = (1.0/(mTimeStep*mTimeStep))*mMassMatrix+mL;
	if(mIntegrationMethod == PROJECTIVE_DYNAMICS)
		FactorizeLDLT(H2ML,mSolver);
	else if(mIntegrationMethod == PROJECTIVE_QUASI_STATIC)
		FactorizeLDLT(mL,mSolver);
}
void
World::
InversionFree(Eigen::VectorXd& x)
{
	//TODO
}
void
World::
IntegratePositionsAndVelocities(const Eigen::VectorXd& x_next)
{
	mVelocities = (1.0/mTimeStep)*(x_next - mPositions);
	mPositions = x_next;
	mVelocities = mDampingCoefficient*mVelocities;
}