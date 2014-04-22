#include "common/common.h"
#include "common/generation.h"
#include "common/parser.h"
#include "common/input_output.h"

bool visualize = 0;
double particle_c = 0; // timestep
bool addParticles = 0;

real gravity = -9.80665;
real timestep = 0.00075;
real seconds_to_simulate = 30;
int num_steps = seconds_to_simulate / timestep;
int max_iter = 100;
real tolerance = 8e-5;

real3 container_size = R3(6.25, 2, 3);
real container_thickness = .2;
real container_height = -1;
real container_friction = 1;

string data_folder = "data/plow";

ChSystemParallelDVI * system_gpu;

ChSharedPtr<ChMaterialSurface> material_shoes, material_chassis;
ChSharedBodyPtr chassis, engine1, engine2, idler;
ChSharedBodyPtr rollers[10];
ChSharedPtr<ChLinkEngine> eng_roller1, eng_roller2;

vector<ChSharedPtr<ChLinkLockRevolute> > revolutes;
vector<ChBody*> shoes;

double idlerPos = 0;
int sim_type = 0;
float mass_multiplier = 10;
float mass_shoe = 1;//2 * mass_multiplier;
float mass_idler = 1;//25 * mass_multiplier;
float mass_chasis = 1;//500 * mass_multiplier;
float mass_sprocket = 1;//25 * mass_multiplier;
float mass_roller = 1;//25 * mass_multiplier;

double L = .09075;
double R = .039;
double H1 = .08;
double H2 = .0739;
double D = .1692;
double W = .015;
int roller_sprocker_counter = 0;
real particle_radius = .04;
real cohesion = 500;

ParticleGenerator<ChSystemParallel>* layer_gen;
ParticleGenerator<ChSystemParallel>* layer_gen_bottom;

ChSharedBodyPtr createTrackShoeM113(ChVector<> position, ChQuaternion<> rotation)
{
	double L = .09075;
	double R = .039;
	double H1 = .08;
	double H2 = .0739;
	double D = .1692;
	double W = .015;

	ChSharedBodyPtr mrigidBody = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

	InitObject(mrigidBody, mass_shoe, position, rotation, material_shoes, true, false, 3, 3);
	AddCollisionGeometry(mrigidBody, SPHERE, 	ChVector<>(R, D, R), 				ChVector<>(L, 0, 0), 					chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	//AddCollisionGeometry(mrigidBody, CYLINDER, 	ChVector<>(R, D, R), 				ChVector<>(L, 0, 0), 					chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	//AddCollisionGeometry(mrigidBody, CYLINDER, 	ChVector<>(R, D, R), 				ChVector<>(-L, 0, 0), 					chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	AddCollisionGeometry(mrigidBody, BOX, 		ChVector<>(L+2*R,(H1+H2)*.5,R*.5), 	ChVector<>(0,H1-(H1+H2)*.5,D+R*.5), 	Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(mrigidBody, BOX, 		ChVector<>(L+2*R,(H1+H2)*.5,R*.5), 	ChVector<>(0,H1-(H1+H2)*.5,-D-R*.5), 	Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(mrigidBody, BOX, 		ChVector<>(L+2*R,W,D), 				ChVector<>(0,-H2+W,0), 					Quaternion(1, 0, 0, 0));
	FinalizeObject(mrigidBody, (ChSystemParallel *) system_gpu);
	//mrigidBody->SetInertiaXX(ChVector<>(.00067, .00082, .00017));

	return mrigidBody;
}

ChSharedBodyPtr createChassisM113(ChVector<> &position)
{
	ChSharedBodyPtr mrigidBody = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

	InitObject(mrigidBody, mass_chasis, position, Quaternion(1, 0, 0, 0), material_chassis, true, false, 4, 4);
	AddCollisionGeometry(mrigidBody, BOX, 		ChVector<>(.1,.1,.1), 	ChVector<>(0,0,0), 	Quaternion(1, 0, 0, 0));
	FinalizeObject(mrigidBody, (ChSystemParallel *) system_gpu);

	return mrigidBody;
}

ChSharedBodyPtr createRollerM113(ChVector<> position, ChSharedBodyPtr mrigidBody1, double springLength, double springK, double springR)
{
	ChSharedBodyPtr mrigidBody = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

	InitObject(mrigidBody, mass_roller, position, chrono::Q_from_AngAxis(0, VECT_X),material_shoes, true, false, 4, 4);
	AddCollisionGeometry(mrigidBody, CYLINDER, ChVector<>(.305,.16,.305), ChVector<>(0, 0, 0), chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	FinalizeObject(mrigidBody, (ChSystemParallel *) system_gpu);

	ChSharedBodyPtr ptr1 = ChSharedBodyPtr(mrigidBody);
	ChSharedBodyPtr ptr2 = ChSharedBodyPtr(mrigidBody1);
	ChSharedPtr<ChLinkLockRevolute> revJoint(new ChLinkLockRevolute);
	revJoint->Initialize(ptr1, ptr2, ChCoordsys<>(position, QUNIT));
	system_gpu->AddLink(revJoint);

	return mrigidBody;
}

ChSharedBodyPtr createSprocketM113(ChVector<> &position, ChSharedBodyPtr mrigidBody1)
{
	ChSharedBodyPtr mrigidBody = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

	InitObject(mrigidBody, mass_sprocket, position, chrono::Q_from_AngAxis(0, VECT_X), material_shoes, true, false, 4, 4);
	AddCollisionGeometry(mrigidBody, CYLINDER, ChVector<>(.254,.16,.254), ChVector<>(0, 0, 0), chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	for (int i = 0; i < 5; i++) {
		AddCollisionGeometry(mrigidBody, BOX, ChVector<>(.28965,.03,.16), ChVector<>(0, 0, 0), chrono::Q_from_AngAxis(i * 36.0 * CH_C_PI / 180.0, VECT_Z));
	}
	FinalizeObject(mrigidBody, (ChSystemParallel *) system_gpu);
	//mrigidBody->SetInertiaXX(ChVector<>(.3717, .3717, .736));

	ChSharedPtr<ChLinkLockRevolute> my_link1(new ChLinkLockRevolute);
	my_link1->Initialize(mrigidBody, mrigidBody1, ChCoordsys<>(position, QUNIT));
	system_gpu->AddLink(my_link1);

	return mrigidBody;
}

ChSharedBodyPtr createIdlerM113(ChVector<> position, ChSharedBodyPtr mrigidBody1, double springLength, double springK, double springR)
{
	ChSharedBodyPtr mrigidBody = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

	InitObject(mrigidBody, mass_idler, position, chrono::Q_from_AngAxis(0, VECT_X), material_shoes, true, false, 4, 4);
	AddCollisionGeometry(mrigidBody, CYLINDER, ChVector<>(.2538,.16,.2538), ChVector<>(0, 0, 0), chrono::Q_from_AngAxis(CH_C_PI / 2, VECT_X));
	FinalizeObject(mrigidBody, (ChSystemParallel *) system_gpu);

	ChSharedPtr<ChLinkLockRevolute> revJoint(new ChLinkLockRevolute);
	revJoint->Initialize(mrigidBody, mrigidBody1, ChCoordsys<>(position, QUNIT));
	system_gpu->AddLink(revJoint);

	return mrigidBody;
}

ChSharedBodyPtr inputTrackModelM113(int nTracks, real3 pinLoc1, ChVector<> &position, ChSharedBodyPtr chassisBody)
{
	ChSharedBodyPtr mrigidBody;
	ChSharedBodyPtr mrigidBodyPrev;
	ChSharedBodyPtr mrigidBodyFirst;

	int bodyNum;
	double pos_x = 0, pos_y = 0, pos_z = 0, rot_x = 0, rot_y = 0, rot_z = 0;

	string temp_data;
	ifstream ifile("pos_M113.dat");
	getline(ifile, temp_data);
	stringstream ss(temp_data);
	for (int i = 0; i < nTracks; i++)
	{
		getline(ifile, temp_data);
		stringstream ss(temp_data);
		ss >> bodyNum >> pos_x >> pos_y >> pos_z >> rot_x >> rot_y >> rot_z;

		mrigidBodyPrev = mrigidBody;
		ChVector<> pos_temp = ChVector<>(pos_x, pos_y, 0) + position;
		ChQuaternion<> temp_Q = chrono::Q_from_AngAxis(rot_z, VECT_Z);
		mrigidBody = createTrackShoeM113(pos_temp, temp_Q);
		//mrigidBody->SetBodyFixed(true);
		shoes.push_back(mrigidBody.get_ptr());

		if (i == 0)
		{
			mrigidBodyFirst = mrigidBody;
		}
		if (i != 0)
		{
			ChCoordsys<> pos1, pos2;
			pos1.pos = Vector(0, 0, 0);
			pos2.pos = Vector(0, 0, 0);
			ChSharedPtr<ChLinkLockRevolute> revJoint(new ChLinkLockRevolute);
			revJoint->Initialize(mrigidBody, mrigidBodyPrev, ChCoordsys<>(mrigidBody->Point_Body2World(ChVector<>(pinLoc1.x, pinLoc1.y, pinLoc1.z)), QUNIT));
			system_gpu->AddLink(revJoint);
			revolutes.push_back(revJoint);
		}
		if (i == (nTracks - 1))
		{
			ChCoordsys<> pos1, pos2;
			pos1.pos = Vector(0, 0, 0);
			pos2.pos = Vector(0, 0, 0);

			ChSharedPtr<ChLinkLockRevolute> revJoint(new ChLinkLockRevolute);
			revJoint->Initialize(mrigidBody, mrigidBodyFirst, ChCoordsys<>(mrigidBodyFirst->Point_Body2World(ChVector<>(pinLoc1.x, pinLoc1.y, pinLoc1.z)), QUNIT));
			system_gpu->AddLink(revJoint);
			revolutes.push_back(revJoint);
		}

	}

	for (int i = 0; i < 5; i++) {
		getline(ifile, temp_data);
		std::stringstream ss(temp_data);
		ss >> bodyNum >> pos_x >> pos_y >> pos_z >> rot_x >> rot_y >> rot_z;
		ChVector<> pos_temp = ChVector<>(pos_x, pos_y, 0) + position;
		rollers[roller_sprocker_counter] = createRollerM113(pos_temp, chassisBody, 0, 0, 0);
		ChQuaternion<> temp;
		temp.Q_from_NasaAngles(ChVector<>(rot_x, rot_y, rot_z));
		rollers[roller_sprocker_counter]->SetRot(temp);
		roller_sprocker_counter++;
	}

	ChSharedBodyPtr my_link1;
	for (int i = 0; i < 1; i++)
	{
		getline(ifile, temp_data);
		std::stringstream ss(temp_data);
		ss >> bodyNum >> pos_x >> pos_y >> pos_z >> rot_x >> rot_y >> rot_z;
		ChVector<> pos_temp = ChVector<>(pos_x, pos_y, 0) + position;
		my_link1 = createSprocketM113(pos_temp, chassisBody);
		ChQuaternion<> temp;
		temp.Q_from_NasaAngles(ChVector<>(rot_x, rot_y, rot_z));
		my_link1->SetRot(temp);
	}
	for (int i = 0; i < 1; i++)
	{
		getline(ifile, temp_data);
		std::stringstream ss(temp_data);
		ss >> bodyNum >> pos_x >> pos_y >> pos_z >> rot_x >> rot_y >> rot_z;
		ChVector<> pos_temp = ChVector<>(pos_x, pos_y, 0) + position;
		ChSharedBodyPtr idler = createIdlerM113(pos_temp, chassisBody, 10, 500, 100);
	}

	return my_link1;
}

template<class T>
void RunTimeStep(T* mSys, const int frame)
{
	engine1->Empty_forces_accumulators();
	engine2->Empty_forces_accumulators();
	if (frame * timestep > .2)
	{
		engine1->SetWvel_loc(ChVector<>(0, 0, 4));
		engine2->SetWvel_loc(ChVector<>(0, 0, 4));
//		for (int i = 0; i < 6; i++)
//		{
//			rollers[i]->SetWvel_loc(ChVector<>(0, 0, 4));
//		}
	}

	((ChSystemParallel*) mSys)->SetAABB(R3(-6.25, -2+container_height, -3), R3(6.25, 2+container_height, 3));
}

int main(int argc, char* argv[])
{
	if (argc > 1) {
		data_folder = argv[1];
		visualize = atoi(argv[2]);
		particle_c = atof(argv[3]); // timestep
		addParticles = atoi(argv[4]);
	}

	// Set up the solver
	system_gpu = new ChSystemParallelDVI;
	system_gpu->SetIntegrationType(ChSystem::INT_ANITESCU);

	((ChLcpSolverParallelDVI*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationNormal(40);
	((ChLcpSolverParallelDVI*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSliding(40);
	((ChLcpSolverParallelDVI*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationSpinning(0);
	((ChLcpSolverParallelDVI*) (system_gpu->GetLcpSolverSpeed()))->SetMaxIterationBilateral(40);
	system_gpu->SetTol(tolerance);
	system_gpu->SetTolSpeeds(tolerance);
	system_gpu->SetMaxPenetrationRecoverySpeed(100);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetTolerance(tolerance);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetCompliance(0);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetContactRecoverySpeed(5);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetSolverType(APGDRS);
	((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->SetWarmStart(false);
	//((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->DoCollision(false);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->SetCollisionEnvelope(particle_radius * .5 * .05);
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBinsPerAxis(I3(60, 20, 30));
	((ChCollisionSystemParallel *) (system_gpu->GetCollisionSystem()))->setBodyPerBin(200, 100);
	system_gpu->Set_G_acc(ChVector<>(0, gravity, 0));
	system_gpu->SetStep(timestep);
	//int threads = 16;
	//system_gpu->SetParallelThreadNumber(threads);
	//omp_set_num_threads(threads);
	// End set up solver

	// Create box to keep everything inside
	ChSharedBodyPtr L = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr R = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr F = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr B = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));
	ChSharedBodyPtr Bottom = ChSharedBodyPtr(new ChBody(new ChCollisionModelParallel));

	ChSharedPtr<ChMaterialSurface> material;
	material = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material->SetFriction(container_friction);

	material->SetCompliance(0);
	material->SetCohesion(0);

	InitObject(L, 100000, Vector(-container_size.x + container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(R, 100000, Vector(container_size.x - container_thickness, container_height - container_thickness, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(F, 100000, Vector(0, container_height - container_thickness, -container_size.z + container_thickness), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(B, 100000, Vector(0, container_height - container_thickness, container_size.z - container_thickness), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);
	InitObject(Bottom, 100000, Vector(0, container_height - container_size.y, 0), Quaternion(1, 0, 0, 0), material, true, true, -20, -20);

	AddCollisionGeometry(L, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(R, BOX, Vector(container_thickness, container_size.y, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(F, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(B, BOX, Vector(container_size.x, container_size.y, container_thickness), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));
	AddCollisionGeometry(Bottom, BOX, Vector(container_size.x, container_thickness, container_size.z), Vector(0, 0, 0), Quaternion(1, 0, 0, 0));

	FinalizeObject(L, (ChSystemParallel *) system_gpu);
	FinalizeObject(R, (ChSystemParallel *) system_gpu);
	FinalizeObject(F, (ChSystemParallel *) system_gpu);
	FinalizeObject(B, (ChSystemParallel *) system_gpu);
	FinalizeObject(Bottom, (ChSystemParallel *) system_gpu);
	// End create box to keep everything inside

	material_shoes = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material_shoes->SetFriction(1);
	material_shoes->SetCompliance(0);
	material_shoes->SetCohesion(0);//-1000);

	material_chassis = ChSharedPtr<ChMaterialSurface>(new ChMaterialSurface);
	material_chassis->SetFriction(1);
	material_chassis->SetCohesion(0);//-1000);

	real height = -1.8;
	real x_offset = 1.25;
	ChVector<double> tankPos = ChVector<>(0,-1.7,0);
	ChVector<> temporary1 = ChVector<>(2.358,0,0)+tankPos;
	ChVector<> temporary2 = ChVector<>(0,0,1.052+.08)+tankPos;
	ChVector<> temporary3 = ChVector<>(0,0,-1.052-.08)+tankPos;

	chassis = createChassisM113(temporary1);
	engine1 = inputTrackModelM113(61, R3(-.09075, 0, 0), temporary2, chassis);
	engine2 = inputTrackModelM113(61, R3(-.09075, 0, 0), temporary3, chassis);

	int3 num_per_dir;
	//num_per_dir = I3(200, 1, 80);
	//num_per_dir = I3(150, 1, 50);
	//num_per_dir = I3(150, 8, 50);

	if(addParticles)
	{
		layer_gen = new ParticleGenerator<ChSystemParallel>((ChSystemParallel *) system_gpu);
		layer_gen->SetDensity(10 * 200);
		layer_gen->SetRadius(R3(particle_radius, particle_radius * .5, particle_radius) * .5);
		layer_gen->material->SetFriction(1);
		layer_gen->material->SetCohesion(particle_c);
		layer_gen->material->SetRollingFriction(0);
		layer_gen->material->SetSpinningFriction(0);
		layer_gen->material->SetCompliance(0);
		layer_gen->AddMixtureType(MIX_ELLIPSOID);

		num_per_dir = I3(100, 16, 100);
		layer_gen->addPerturbedVolumeMixture(R3(-2.5, -1.9, 0), I3(num_per_dir.x, num_per_dir.y, num_per_dir.z), R3(.01, .01, .01), R3(0, 0, 0));

		layer_gen_bottom = new ParticleGenerator<ChSystemParallel>((ChSystemParallel *) system_gpu);
		layer_gen_bottom->SetDensity(10 * 200);
		layer_gen_bottom->SetRadius(R3(particle_radius, particle_radius * .5, particle_radius));
		layer_gen_bottom->material->SetFriction(1);
		layer_gen_bottom->material->SetCohesion(500 * particle_c);
		layer_gen_bottom->material->SetRollingFriction(0);
		layer_gen_bottom->material->SetSpinningFriction(0);
		layer_gen_bottom->material->SetCompliance(0);
		layer_gen_bottom->AddMixtureType(MIX_ELLIPSOID);
		layer_gen_bottom->addPerturbedVolumeMixture(R3(0, -2.7, 0), I3(130, 3, 50), R3(.01, .01, .01), R3(0, 0, 0));
	}

//=========================================================================================================
//Rendering specific stuff:
	if(visualize)
	{
		ChOpenGLManager * window_manager = new ChOpenGLManager();
		ChOpenGL openGLView(window_manager, system_gpu, 800, 600, 0, 0, "Test_Solvers");
		openGLView.render_camera->camera_position = glm::vec3(0, -5, -10);
		openGLView.render_camera->camera_look_at = glm::vec3(0, -5, 0);
		openGLView.render_camera->camera_scale = .075;
		openGLView.SetCustomCallback(RunTimeStep);
		openGLView.StartSpinning(window_manager);
		window_manager->CallGlutMainLoop();
	}
//=========================================================================================================
	int file = 0;
	for (int i = 0; i < num_steps; i++) {
		system_gpu->DoStepDynamics(timestep);
		double TIME = system_gpu->GetChTime();
		double STEP = system_gpu->GetTimerStep();
		double BROD = system_gpu->GetTimerCollisionBroad();
		double NARR = system_gpu->GetTimerCollisionNarrow();
		double LCP = system_gpu->GetTimerLcp();
		double UPDT = system_gpu->GetTimerUpdate();
		double RESID = ((ChLcpSolverParallelDVI *) (system_gpu->GetLcpSolverSpeed()))->GetResidual();
		int BODS = system_gpu->GetNbodies();
		int CNTC = system_gpu->GetNcontacts();
		int REQ_ITS = ((ChLcpSolverParallelDVI*) (system_gpu->GetLcpSolverSpeed()))->GetTotalIterations();

		printf("%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7.4f|%7d|%7d|%7d|%7.4f\n", TIME, STEP, BROD, NARR, LCP, UPDT, BODS, CNTC, REQ_ITS, RESID);

		int save_every = 1.0 / timestep / 60.0;     //save data every n steps
		if (i % save_every == 0) {
			stringstream ss;
			cout << "Frame: " << file << endl;
			ss << data_folder << "/" << file << ".txt";
			DumpAllObjectsWithGeometryPovray(system_gpu, ss.str());
			file++;
		}
		RunTimeStep(system_gpu, i);
	}
}
