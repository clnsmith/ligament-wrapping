// Copyright (c) 2011, Eftychios Sifakis.
// Distributed under the FreeBSD license (see license.txt)

#include <PhysBAM_Geometry/Geometry_Particles/REGISTER_GEOMETRY_READ_WRITE.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/SEGMENTED_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>


#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>

#include "SIMULATION_LAYOUT.h"
using namespace PhysBAM;

template<class T>
SIMULATION_LAYOUT<T>::SIMULATION_LAYOUT(const STREAM_TYPE stream_type_input)
    :stream_type(stream_type_input),
     n(21),youngs_modulus(1e1),damping_coefficient(2e0),
    wire_mass(1.),wire_restlength(1.),
    collection(particles),
    number_of_frames(100),frame_time(.05),CFL_number(0.03),
    sphere_radius((T).25),sphere_position(TV(.5,.05,0)),
    collision_stiffness(1e3)
{
    Initialize_Geometry_Particle();Initialize_Read_Write_Structures();
}

template<class T>
void SIMULATION_LAYOUT<T>::Initialize()
{
    particles.Store_Velocity();

    // Wire 1
    wire_curve=SEGMENTED_CURVE<TV>::Create(particles);
    wire_curve->mesh.Initialize_Straight_Mesh(n);particles.array_collection->Add_Elements(n);
    for(int p=1;p<=n;p++) particles.X(p)=TV((T)(p-1)/(T)(n-1),0,0);

    wire_particles=FREE_PARTICLES<TV>::Create(particles);
    for(int p=1;p<=n;p++) wire_particles->nodes.Append(p);

    collection.Add_Structure(wire_curve);collection.Add_Structure(wire_particles);

    mass.Resize(n);mass.Fill(wire_mass/(T)n);
    restlength.Resize(n-1);restlength.Fill(wire_restlength/(T)(n-1));

    // Wire 2
    wire_curve2=SEGMENTED_CURVE<TV>::Create(particles);
    wire_curve2->mesh.Initialize_Straight_Mesh(n);particles.array_collection->Add_Elements(n);
    for(int p=1;p<=n;p++) particles.X(p)=TV((T)(p-1)/(T)(n-1),0,0.25);

    wire_particles2=FREE_PARTICLES<TV>::Create(particles);
    for(int p=1;p<=n;p++) wire_particles2->nodes.Append(p);

    collection.Add_Structure(wire_curve2);collection.Add_Structure(wire_particles2);

    mass.Resize(n);mass.Fill(wire_mass/(T)n);
    restlength.Resize(n-1);restlength.Fill(wire_restlength/(T)(n-1));

#if 0
    TETRAHEDRALIZED_VOLUME<T>* ball_volume=TETRAHEDRALIZED_VOLUME<T>::Create();
    FILE_UTILITIES::Read_From_File(stream_type,"sphere.tet",*ball_volume);
    ball_volume->Initialize_Triangulated_Surface();
    ball_volume->triangulated_surface->Discard_Valence_Zero_Particles_And_Renumber();
    FILE_UTILITIES::Write_To_File(stream_type,"sphere.tri",*ball_volume->triangulated_surface);
    exit(0);
#endif

// Sphere

    TRIANGULATED_SURFACE<T>* input_collision_surface=TRIANGULATED_SURFACE<T>::Create();
    FILE_UTILITIES::Read_From_File(stream_type,"sphere.tri",*input_collision_surface);
    for(int p=1;p<=input_collision_surface->particles.array_collection->Size();p++)
        input_collision_surface->particles.X(p)=input_collision_surface->particles.X(p)*sphere_radius+sphere_position;
    TRIANGULATED_SURFACE<T>* collision_surface=dynamic_cast<TRIANGULATED_SURFACE<T>*>(
        input_collision_surface->Append_Particles_And_Create_Copy(particles));
    collection.Add_Structure(collision_surface);

    wire_curve->Update_Number_Nodes();
    wire_particles->Update_Number_Nodes();

    wire_curve2->Update_Number_Nodes();
    wire_particles2->Update_Number_Nodes();

    collision_surface->Update_Number_Nodes();    
    
    // ball_volume=TETRAHEDRALIZED_VOLUME<T>::Create(particles);
    // FILE_UTILITIES::Read_From_File(stream_type,"sphere.tet",*ball_volume);
    // LOG::cout<<"Initializing sphere geometry ..."<<std::endl;
    // LOG::cout<<"  Particles="<<ball_volume->particles.array_collection->Size()<<std::endl;
    // LOG::cout<<"  Tetrahedrons="<<ball_volume->mesh.elements.m<<std::endl;

}

template<class T>
void SIMULATION_LAYOUT<T>::Add_Elastic_Forces(ARRAY_VIEW<const TV> X,ARRAY_VIEW<TV> force) const
{
    //Wire 1 Elastic
    for(int s=1;s<=wire_curve->mesh.elements.m;s++){
        int p1,p2;wire_curve->mesh.elements(s).Get(p1,p2);
        TV X1=X(p1),X2=X(p2);
        TV f=-(youngs_modulus/restlength(s))*(X1-X2);
        force(p1)+=f;force(p2)-=f;
    }

    //Wire 2 Elastic
    for(int s=1;s<=wire_curve2->mesh.elements.m;s++){
        int p1,p2;wire_curve2->mesh.elements(s).Get(p1,p2);
        TV X1=X(p1),X2=X(p2);
        TV f=-(youngs_modulus/restlength(s))*(X1-X2);
        force(p1)+=f;force(p2)-=f;
    }

    //Sphere Contact
    for(int p=1;p<=2*n;p++){
        const TV& X=particles.X(p);
        T depth=(X-sphere_position).Magnitude()-sphere_radius;
        if(depth>0) continue;
        TV normal=(X-sphere_position).Normalized();
        TV f=(-collision_stiffness*depth)*normal;
        force(p)+=f;
    }

}

template<class T>
void SIMULATION_LAYOUT<T>::Add_Damping_Forces(ARRAY_VIEW<const TV> X,ARRAY_VIEW<const TV> V,
    ARRAY_VIEW<TV> force) const
{
    //Wire 1
    for(int s=1;s<=wire_curve->mesh.elements.m;s++){
        int p1,p2;wire_curve->mesh.elements(s).Get(p1,p2);
        TV V1=V(p1),V2=V(p2);
        TV f=-(damping_coefficient/restlength(s))*(V1-V2);
        force(p1)+=f;force(p2)-=f;

    //Wire 2
    for(int s=1;s<=wire_curve2->mesh.elements.m;s++){
        int p1,p2;wire_curve2->mesh.elements(s).Get(p1,p2);
        TV V1=V(p1),V2=V(p2);
        TV f=-(damping_coefficient/restlength(s))*(V1-V2);
        force(p1)+=f;force(p2)-=f;
    }
}

template<class T>
void SIMULATION_LAYOUT<T>::Add_Force_Differential(ARRAY_VIEW<const TV> X,ARRAY_VIEW<const TV> dX,
    ARRAY_VIEW<TV> dforce) const
{
    //Wire 1
    for(int s=1;s<=wire_curve->mesh.elements.m;s++){
        int p1,p2;wire_curve->mesh.elements(s).Get(p1,p2);
        TV dX1=dX(p1),dX2=dX(p2);
        TV df=-(youngs_modulus/restlength(s))*(dX1-dX2);
        dforce(p1)+=df;dforce(p2)-=df;
    }

    //Wire 2
    for(int s=1;s<=wire_curve2->mesh.elements.m;s++){
        int p1,p2;wire_curve2->mesh.elements(s).Get(p1,p2);
        TV dX1=dX(p1),dX2=dX(p2);
        TV df=-(youngs_modulus/restlength(s))*(dX1-dX2);
        dforce(p1)+=df;dforce(p2)-=df;
    }
}

template<class T>
void SIMULATION_LAYOUT<T>::Add_External_Forces(ARRAY_VIEW<TV> force) const
{
    //Gravity on Wire
//    for(int p=1;p<=n;p++)
//        force(p)-=TV::Axis_Vector(2)*mass(p)*9.81;
}

template<class T>
T SIMULATION_LAYOUT<T>::Maximum_Dt() const
{
    return frame_time*CFL_number;
}

template<class T>
void SIMULATION_LAYOUT<T>::Write_Output(const int frame) const
{
    FILE_UTILITIES::Create_Directory("output/"+STRING_UTILITIES::Value_To_String(frame));
    collection.Write(stream_type,"output",frame,0,true);
}

template<class T>
void SIMULATION_LAYOUT<T>::Set_Kinematic_Positions(const T time,ARRAY_VIEW<TV> X) const
{
    // Fix the 2 Wire end points to remain immovable
    //Wire 1
    X(1)=TV();
    X(n)=TV::Axis_Vector(1);

    //Wire 2
    //X(n+1)=TV(0,0,0.25);
    //X(2*n)=TV(1,0,0.25);
    X(n+1)=TV();
    X(2*n)=TV(0,0,1);
}

template<class T>
void SIMULATION_LAYOUT<T>::Add_Kinematic_Positions(const T time,const T factor,ARRAY_VIEW<TV> X) const
{
    // Fix the 2 extremes to remain immovable
    
    //Wire 1
    X(1)+=TV();
    X(n)+=TV::Axis_Vector(1)*factor;

    //Wire 2
    //X(n+1)=TV(0,0,0.25);
    //X(2*n)=TV(1,0,0.25);

    X(n+1)=TV();
    X(2*n)=TV(0,0,1)*factor;
}

template<class T>
void SIMULATION_LAYOUT<T>::Set_Kinematic_Velocities(const T time,ARRAY_VIEW<TV> V) const
{
    // Fix the 2 extremes to remain immovable
    //Wire1
    V(1)=TV();
    V(n)=TV();

    //Wire2
    V(n+1)=TV();
    V(2*n)=TV();

    //Sphere
    for(int p=2*n+1;p<=particles.array_collection->Size();p++)
        V(p)=TV();
}

template<class T>
void SIMULATION_LAYOUT<T>::Clear_Values_Of_Kinematic_Particles(ARRAY_VIEW<TV> array) const
{
    // Fix the 2 extremes to remain immovable
    //Wire 1
    array(1)=TV();
    array(n)=TV();

    //Wire 2
    array(n+1)=TV();
    array(2*n)=TV();

    //Sphere
    for(int p=2*n+1;p<=particles.array_collection->Size();p++)
        array(p)=TV();
}

template class SIMULATION_LAYOUT<float>;
