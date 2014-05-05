// Copyright (c) 2011, Eftychios Sifakis.
// Distributed under the FreeBSD license (see license.txt)

#pragma once

namespace PhysBAM{

template<class TV> class SEGMENTED_CURVE;
template<class TV> class FREE_PARTICLES;

template<class T>
class SIMULATION_LAYOUT
{
public:
    typedef VECTOR<T,3> TV;

    const STREAM_TYPE stream_type;

    const int n;                 // Number of particles in wire mesh
    const T youngs_modulus;      // Elasticity and damping coefficients
    const T damping_coefficient;
    const T wire_mass;           // Mass and length for entire wire
    const T wire_restlength;
    ARRAY<T> mass;               // Mass (per each particle)
    ARRAY<T> restlength;         // Restlength (per each spring)
    
    GEOMETRY_PARTICLES<TV> particles;
    DEFORMABLE_GEOMETRY_COLLECTION<TV> collection;

    const int number_of_frames;  // Total number of frames
    const T frame_time;          // Frame (snapshot) interval
    const T CFL_number;          // CFL number (not to exceed 1)

    SEGMENTED_CURVE<TV>* wire_curve;
    FREE_PARTICLES<TV>* wire_particles;

    SEGMENTED_CURVE<TV>* wire_curve2;
    FREE_PARTICLES<TV>* wire_particles2;

    const T sphere_radius;
    const TV sphere_position;

    const T collision_stiffness;
    
    SIMULATION_LAYOUT(const STREAM_TYPE stream_type_input);
    void Initialize();
    void Add_Elastic_Forces(ARRAY_VIEW<const TV> X,ARRAY_VIEW<TV> force) const;
    void Add_Damping_Forces(ARRAY_VIEW<const TV> X,ARRAY_VIEW<const TV> V,
        ARRAY_VIEW<TV> force) const;
    void Add_Force_Differential(ARRAY_VIEW<const TV> X,ARRAY_VIEW<const TV> dX,
        ARRAY_VIEW<TV> dforce) const;
    void Add_External_Forces(ARRAY_VIEW<TV> force) const;
    T Maximum_Dt() const;
    void Write_Output(const int frame) const;
    void Set_Kinematic_Positions(const T time,ARRAY_VIEW<TV> X) const;
    void Add_Kinematic_Positions(const T time,const T factor,ARRAY_VIEW<TV> X) const;
    void Set_Kinematic_Velocities(const T time,ARRAY_VIEW<TV> V) const;
    void Clear_Values_Of_Kinematic_Particles(ARRAY_VIEW<TV> array) const;
};

}
