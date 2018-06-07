/**
 * @file scene.hpp
 * @brief Class definitions for scenes.
 *
 */

#ifndef _462_SCENE_SCENE_HPP_
#define _462_SCENE_SCENE_HPP_

#include "math/vector.hpp"
#include "math/quaternion.hpp"
#include "math/matrix.hpp"
#include "math/camera.hpp"
#include "scene/material.hpp"
#include "scene/mesh.hpp"
//#include "scene/bvh.hpp"
#include "ray.hpp"
#include <string>
#include <vector>
#include <cfloat>

namespace _462 {

// lower time bound 
const real_t t0 = 0.00015462;
// upper time bound
const real_t t1 = 15462999.0;
const real_t neg_infinity = -999999.9;
const real_t pos_infinity =  999999.9;

/**
 * Details of the neighbor photon stored in max heap
 */
struct Neighbor
{
    Neighbor() {
        sq_dis = -1.0;
        i = 0;
    };
    Neighbor(real_t s, unsigned e) {
        sq_dis = s;
        i = e;
    };
    // squared distance from photon to intersection
    real_t sq_dis;
    // index of photon in map
    unsigned i;
};

/**
 * Details of the photon emitted from light source
 */
struct Photon
{
    Photon() {};
    Photon(char f) {
        flag = f;
    };
    Photon(Color3 i){ 
        c = i;
    };
    // indicate the splitting axis or leaf node
    char flag;
    // position of photon
    Vector3 p;
    // incident direction of photon ray
    Vector3 dir;
    // intensity of photon
    Color3 c;
};

/**
 * Details of the surface closest to the eye ray
 */
struct Intersection
{
    // Iniitialize variables
    Intersection() {
        local_t = t1;
        material.ambient = Color3::Black();
        material.diffuse = Color3::Black();
        material.specular = Color3::Black();
        material.refractive_index = 0.0;
        texture_color = Color3::Black();
        world_intersection = Vector3::Zero();
        world_normal = Vector3::Zero();
    };

    real_t local_t;
    Material material;
    Color3 texture_color;
    // location of the intersection in world space
    Vector3 world_intersection;
    // normal vector on the surface of intersected object
    Vector3 world_normal;
};

/**
 * Details of the bounding box of a geometry
 */
struct Bound
{
    Bound() {
        min_x = pos_infinity;
        max_x = neg_infinity;
        min_y = pos_infinity;
        max_y = neg_infinity;
        min_z = pos_infinity;
        max_z = neg_infinity;
    };
    // the minimum and maximum of the coordinate for the geometry 
    // so that the box completely contains the geometry
    real_t min_x, max_x, min_y, max_y, min_z, max_z;
};   

class Geometry
{
public:
    Geometry();
    virtual ~Geometry();

    /*
       World transformation are applied in the following order:
       1. Scale
       2. Orientation
       3. Position
    */

    // The world position of the object.
    Vector3 position;

    // The world orientation of the object.
    // Use Quaternion::to_matrix to get the rotation matrix.
    Quaternion orientation;

    // The world scale of the object.
    Vector3 scale;

    // Transformation matrix
    Matrix4 mat;

    // Inverse transformation matrix
    Matrix4 invMat;
    // Normal transformation matrix
    Matrix3 normMat;

    // Bounding box of the geometry 
    Bound bound;
     
    // Transfrom the ray from world space to local space
    Ray world_to_local_ray(const Ray& r) const;
    
    // Calculate the location of intersection in world space
    Vector3 local_to_world_intersection(const Ray& local_r, 
                                        real_t local_t) const;
    
    // Calculate the normal in world space
    Vector3 local_to_world_normal(Vector3 local_n) const;
    
    //Calculate the texture color given the tecture coordinate
    Color3 texture_lookup(Vector2 tex_coord, const Material *material) const;

    /**
     * Renders this geometry using OpenGL in the local coordinate space.
     */
    virtual void render() const = 0;
    
    /**
     * Checks whether the given ray hits this geometry in local space.
     */    
    virtual void intersection_test(const Ray& r, Intersection& inter, 
                                   real_t *closest_found_t) const = 0; 

    /**
     * Checks whether the given ray hits this geometry in local space.
     */    
    virtual void init_bound(Bound& bound) const = 0; 

    /**
     * Initialize matrices and bound.
     */    
    bool initialize();

    /**
     * Checks whether the given ray hits this geometry's bound.
     */    
    bool intersection_bound_test(const Ray& r, real_t *closest_found_t) const;
};


struct SphereLight
{
    struct Attenuation
    {
        real_t constant;
        real_t linear;
        real_t quadratic;
    };

    SphereLight();

    bool intersect(const Ray& r, real_t& t);

    // The position of the light, relative to world origin.
    Vector3 position;
    // The color of the light (both diffuse and specular)
    Color3 color;
    // attenuation
    Attenuation attenuation;
    real_t radius;
};

/**
 * The container class for information used to render a scene composed of
 * Geometries.
 */
class Scene
{
public:

    /// the camera
    Camera camera;
    /// allow the user to specify the side length of the camera len
    /// to create the depth of field effect for ray tracing
    int camera_len_length;
    /// the background color
    Color3 background_color;
    /// the amibient light of the scene
    Color3 ambient_light;
    /// the refraction index of air
    real_t refractive_index;
    /// allow the user to specify whether light should attenuate in
    /// dielectrics according to Beer's Law
    bool Beers_Law;

    /// Creates a new empty scene.
    Scene();

    /// Destroys this scene. Invokes delete on everything in geometries.
    ~Scene();

    bool initialize();

    // accessor functions
    Geometry* const* get_geometries() const;
    size_t num_geometries() const;
    const SphereLight* get_lights() const;
    size_t num_lights() const;
    Material* const* get_materials() const;
    size_t num_materials() const;
    Mesh* const* get_meshes() const;
    size_t num_meshes() const;

    /// Clears the scene, and invokes delete on everything in geometries.
    void reset();

    // functions to add things to the scene
    // all pointers are deleted by the scene upon scene deconstruction.
    void add_geometry( Geometry* g );
    void add_material( Material* m );
    void add_mesh( Mesh* m );
    void add_light( const SphereLight& l );

private:

    typedef std::vector< SphereLight > SphereLightList;
    typedef std::vector< Material* > MaterialList;
    typedef std::vector< Mesh* > MeshList;
    typedef std::vector< Geometry* > GeometryList;

    // list of all lights in the scene
    SphereLightList point_lights;
    // all materials used by geometries
    MaterialList materials;
    // all meshes used by models
    MeshList meshes;
    // list of all geometries. deleted in dctor, so should be allocated on heap.
    GeometryList geometries;

private:

    // no meaningful assignment or copy
    Scene(const Scene&);
    Scene& operator=(const Scene&);

};

} /* _462 */

#endif /* _462_SCENE_SCENE_HPP_ */
