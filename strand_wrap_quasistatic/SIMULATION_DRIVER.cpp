// Copyright (c) 2011, Eftychios Sifakis.
// Distributed under the FreeBSD license (see license.txt)

#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>

#include "SIMULATION_DRIVER.h"
#include "CG_VECTOR.h"
#include "CG_SYSTEM.h"
using namespace PhysBAM;

template<class T>
void SIMULATION_DRIVER<T>::Simulate_Frame(const int frame)
{
    T frame_end_time=layout.frame_time*(T)frame,dt_max,dt;

    for(;time<frame_end_time;time+=dt){
        dt_max=layout.Maximum_Dt();
        dt=std::min(dt_max,(T)1.001*(frame_end_time-time));
        Simulate_Time_Step(time,dt);}
}

template<class T>
void SIMULATION_DRIVER<T>::Simulate_Time_Step(const T time,const T dt)
{
    layout.Set_Kinematic_Positions(time,layout.particles.X);
    layout.Set_Kinematic_Velocities(time,layout.particles.V);

    // Construct right-hand-side
    const int number_of_particles=layout.particles.array_collection->Size();
    ARRAY<TV> rhs(number_of_particles);
    layout.Add_Elastic_Forces(layout.particles.X,rhs);
    layout.Add_External_Forces(rhs);
    for(int p=1;p<=number_of_particles;p++)
        rhs(p)+=(layout.mass(p)/dt)*layout.particles.V(p);

    // Initial guess - zero position change
    ARRAY<TV> dX(number_of_particles);

    // Temporary vectors, required by Conjugate Gradients
    ARRAY<TV> temp_q(number_of_particles),temp_s(number_of_particles),temp_r(number_of_particles),
        temp_k(number_of_particles),temp_z(number_of_particles);

    // Encapsulate all vectors in CG-mandated format
    CG_VECTOR<T> cg_x(dX),cg_b(rhs),
        cg_q(temp_q),cg_s(temp_s),cg_r(temp_r),cg_k(temp_k),cg_z(temp_z);

    // Generate CG-formatted system object
    CG_SYSTEM<T> cg_system(layout,time,dt);

    // Generate Conjugate Gradients solver object
    CONJUGATE_GRADIENT<T> cg;
    cg.print_residuals=true;
    cg.print_diagnostics=true;

    // Solve linear system using CG
    cg.Solve(cg_system,
        cg_x,cg_b,cg_q,cg_s,cg_r,cg_k,cg_z,
        1e-4,0,50);

    layout.Clear_Values_Of_Kinematic_Particles(dX);

    layout.particles.X+=dX;        // Update particle positions and velocities
    dX*=(1./dt);
    layout.particles.V=dX; // Update particle positions and velocities

    layout.Set_Kinematic_Positions(time+dt,layout.particles.X);
    layout.Set_Kinematic_Velocities(time+dt,layout.particles.V);
}

template class SIMULATION_DRIVER<float>;

