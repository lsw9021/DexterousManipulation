#include "Tensor3333.h"
#include "World.h"
#include "Constraint/Constraint.h"
#include "Mesh/OBJMesh.h"
#include "Mesh/RectangularMesh.h"
#include "Mesh/DiamondMesh.h"

#define MAKE_SHARED_WEAK_TYPEDEF( X )\
	class X ;\
	typedef std::shared_ptr< X >       X ## Ptr;\
	typedef std::shared_ptr< const X > Const ## X ## Ptr;\
	typedef std::weak_ptr< X >         Weak ## X ## Ptr;\
	typedef std::weak_ptr< const X >   WeakConst ## X ## Ptr;

namespace FEM
{
	//Constraint
	MAKE_SHARED_WEAK_TYPEDEF(Cst);
	MAKE_SHARED_WEAK_TYPEDEF(AttachmentCst);
	MAKE_SHARED_WEAK_TYPEDEF(CorotateFEMCst);
	MAKE_SHARED_WEAK_TYPEDEF(LinearMuscleCst);
	
	//Mesh	
	MAKE_SHARED_WEAK_TYPEDEF(Mesh);
	MAKE_SHARED_WEAK_TYPEDEF(RectangularMesh);
	MAKE_SHARED_WEAK_TYPEDEF(DiamondMesh);

	//World
	MAKE_SHARED_WEAK_TYPEDEF(World);
};