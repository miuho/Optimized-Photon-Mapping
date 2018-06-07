/**
 * @file scene.cpp
 * @brief Function definitions for scenes.
 *
 * @author Eric Butler (edbutler)
 * @author Kristin Siu (kasiu)
 */

#include "scene/scene.hpp"

namespace _462 {

Geometry::Geometry():
    position(Vector3::Zero()),
    orientation(Quaternion::Identity()),
    scale(Vector3::Ones())
{

}

Geometry::~Geometry() { }

bool Geometry::initialize()
{
	make_inverse_transformation_matrix(&invMat, position, orientation, scale);
	make_transformation_matrix(&mat, position, orientation, scale);
	make_normal_matrix(&normMat, mat);
	return true;
}

Ray Geometry::world_to_local_ray(const Ray& r) const {
    // transform the eye position and direction vector with
    // the inverse transformation matrix
    return Ray(invMat.transform_point(r.e), 
               invMat.transform_vector(r.d));
}

Vector3 Geometry::local_to_world_intersection(const Ray& local_r,
                                              real_t local_t) const {
    // apply M on the local intersection point to get the world
    // intersection point
    return mat.transform_point(local_r.e + local_t * local_r.d);
}

Vector3 Geometry::local_to_world_normal(Vector3 local_n) const {
    // multiply the local normal vector with normal matrix
    return normalize(normMat * local_n);
}

Color3 Geometry::texture_lookup(Vector2 tex_coord, 
                                const Material* material) const {
    int width, height;
    material->get_texture_size(&width, &height);
    // check if 0 to avoid floating point exception in modulo
    if (width == 0 || height == 0) { 
        // return white if there is no texture
        return Color3::White();
    }
    else {
        // modulo by width and height to avoid seg fault in accessing array
        int x = (int)(width * tex_coord.x) % width;
        int y = (int)(height * tex_coord.y) % height;
        return material->get_texture_pixel(x, y);
    }
}

inline real_t min2(real_t a, real_t b) {
    return (a < b) ? a : b;
}

inline real_t min4(real_t a, real_t b, real_t c, real_t d) {
    return min2(min2(min2(a, b), c), d);
}

inline real_t max2(real_t a, real_t b) {
    return (a > b) ? a : b;
}

inline real_t max4(real_t a, real_t b, real_t c, real_t d) {
    return max2(max2(max2(a, b), c), d);
}

bool Geometry::intersection_bound_test(const Ray& r, real_t *closest_found_t) const {
    real_t t_min, t_max, t_xmin, t_xmax, t_ymin, t_ymax, t_zmin, t_zmax;
    real_t x_a = 1.0/r.d.x, y_a = 1.0/r.d.y, z_a = 1.0/r.d.z;
    real_t x_e = r.e.x, y_e = r.e.y, z_e = r.e.z;
    
    // calculate t interval in x-axis
    if (x_a >= 0) {
        t_xmin = (bound.min_x - x_e) * x_a;
        t_xmax = (bound.max_x - x_e) * x_a;
    }
    else {
        t_xmin = (bound.max_x - x_e) * x_a;
        t_xmax = (bound.min_x - x_e) * x_a;
    }

    // calculate t interval in y-axis
    if (y_a >= 0) {
        t_ymin = (bound.min_y - y_e) * y_a;
        t_ymax = (bound.max_y - y_e) * y_a;
    }
    else {
        t_ymin = (bound.max_y - y_e) * y_a;
        t_ymax = (bound.min_y - y_e) * y_a;
    }

    // calculate t interval in z-axis
    if (z_a >= 0) {
        t_zmin = (bound.min_z - z_e) * z_a;
        t_zmax = (bound.max_z - z_e) * z_a;
    }
    else {
        t_zmin = (bound.max_z - z_e) * z_a;
        t_zmax = (bound.min_z - z_e) * z_a;
    }

    // find if there an intersection among three t intervals
    t_min = max4(t_xmin, t_ymin, t_zmin, t0);
    t_max = min4(t_xmax, t_ymax, t_zmax, *closest_found_t);
    
    return (t_min <= t_max);
}

SphereLight::SphereLight():
    position(Vector3::Zero()),
    color(Color3::White()),
	radius(real_t(0))
{
    attenuation.constant = 1;
    attenuation.linear = 0;
    attenuation.quadratic = 0;
}

Scene::Scene()
{
    reset();
}

Scene::~Scene()
{
    reset();
}

bool Scene::initialize()
{
	bool res = true;
	for (unsigned int i = 0; i < num_geometries(); i++) {
		res &= geometries[i]->initialize();
                geometries[i]->bound = Bound();
                geometries[i]->init_bound(geometries[i]->bound);
        }
	return res;
}

Geometry* const* Scene::get_geometries() const
{
    return geometries.empty() ? NULL : &geometries[0];
}

size_t Scene::num_geometries() const
{
    return geometries.size();
}

const SphereLight* Scene::get_lights() const
{
    return point_lights.empty() ? NULL : &point_lights[0];
}

size_t Scene::num_lights() const
{
    return point_lights.size();
}

Material* const* Scene::get_materials() const
{
    return materials.empty() ? NULL : &materials[0];
}

size_t Scene::num_materials() const
{
    return materials.size();
}

Mesh* const* Scene::get_meshes() const
{
    return meshes.empty() ? NULL : &meshes[0];
}

size_t Scene::num_meshes() const
{
    return meshes.size();
}

void Scene::reset()
{
    for ( GeometryList::iterator i = geometries.begin(); i != geometries.end(); ++i ) {
        delete *i;
    }
    for ( MaterialList::iterator i = materials.begin(); i != materials.end(); ++i ) {
        delete *i;
    }
    for ( MeshList::iterator i = meshes.begin(); i != meshes.end(); ++i ) {
        delete *i;
    }

    geometries.clear();
    materials.clear();
    meshes.clear();
    point_lights.clear();

    camera = Camera();

    background_color = Color3::Black();
    ambient_light = Color3::Black();
    refractive_index = 1.0;
}

void Scene::add_geometry( Geometry* g )
{
    geometries.push_back( g );
}

void Scene::add_material( Material* m )
{
    materials.push_back( m );
}

void Scene::add_mesh( Mesh* m )
{
    meshes.push_back( m );
}

void Scene::add_light( const SphereLight& l )
{
    point_lights.push_back( l );
}


} /* _462 */

