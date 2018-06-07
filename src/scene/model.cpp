/**
 * @file model.cpp
 * @brief Model class
 *
 * @author Eric Butler (edbutler)
 * @author Zeyang Li (zeyangl)
 */

#include "scene/model.hpp"
#include "scene/material.hpp"
#include "application/opengl.hpp"
#include "scene/triangle.hpp"
#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <sstream>


namespace _462 {

Model::Model() : mesh( 0 ), material( 0 ) { }
Model::~Model() { }

void Model::init_bound(Bound& bound) const {
    size_t z;
    // find the minimum and maximum coordinate of all trianges in the model
    for (z = 0; z < mesh->num_vertices(); z++) {
        Vector3 p = mat.transform_point((mesh->get_vertices())[z].position);
        if (p.x < bound.min_x)
            bound.min_x = p.x;
        if (p.x > bound.max_x)
            bound.max_x = p.x;        
        if (p.y < bound.min_y)
            bound.min_y = p.y;
        if (p.y > bound.max_y)
            bound.max_y = p.y;      
        if (p.z < bound.min_z)
            bound.min_z = p.z;
        if (p.z > bound.max_z)
            bound.max_z = p.z;      
    }
} 

void Model::intersection_test(const Ray& r, Intersection& inter,
                              real_t *closest_found_t) const {
    // checking the bound first 
    if (intersection_bound_test(r, closest_found_t) == false)
        return;

    // perform the intersection test in object's local space rather than
    // world space, and so the center of the object is <0, 0, 0>
    Ray local_ray = world_to_local_ray(r);
    Vector3 local_e = local_ray.e;
    Vector3 local_d = local_ray.d;   
    
    size_t z, hitobject;
    unsigned int i0, i1, i2;
    real_t a, b, c, d, e, f, g, h, i, j, k, l, M, t, min_t, 
           Gamma, min_Gamma, Beta, min_Beta, min_Alpha;
    MeshVertex v0, v1, v2, min_v0, min_v1, min_v2;
    // set the upper bound for local t
    min_t = *closest_found_t;
    g = local_d.x;
    h = local_d.y;
    i = local_d.z;
    for (z = 0; z < mesh->num_triangles(); z++) {
        i0 = (mesh->get_triangles())[z].vertices[0];
        i1 = (mesh->get_triangles())[z].vertices[1];
        i2 = (mesh->get_triangles())[z].vertices[2];
        v0 = (mesh->get_vertices())[i0];
        v1 = (mesh->get_vertices())[i1];
        v2 = (mesh->get_vertices())[i2];
    
         // Using Cramer's rule to solve for Beta, Gamma, and t
        a = v0.position.x - v1.position.x;
        b = v0.position.y - v1.position.y;
        c = v0.position.z - v1.position.z;
        d = v0.position.x - v2.position.x;
        e = v0.position.y - v2.position.y;
        f = v0.position.z - v2.position.z;
        j = v0.position.x - local_e.x;
        k = v0.position.y - local_e.y;
        l = v0.position.z - local_e.z;
        M = a * (e * i - h * f) + b * (g * f - d * i) + 
            c * (d * h - e * g);
    
        t = -(f * (a * k - j * b) + e *(j * c - a * l) + 
            d * (b * l - k * c))/M;
        if (t < t0 || t > min_t) {
            continue;
        } 
    
        Gamma = (i * (a * k - j * b) + h * (j * c - a * l) + 
                g * (b * l - k * c))/M;
        if (Gamma < 0.0 || Gamma > 1.0) {
            continue;
        }

        Beta = (j * (e * i - h * f) + k * (g * f - d * i) + 
               l * (d * h - e * g))/M;
        if (Beta < 0.0 || Beta > (1.0 - Gamma)) {
            continue;
        }
        
        hitobject = z;
        min_v0 = v0;
        min_v1 = v1;
        min_v2 = v2;
        min_t = t;
        min_Gamma = Gamma;
        min_Beta = Beta;
    }
    // check if we have a closer t
    if (min_t == *closest_found_t) {
        return;
    }

    min_Alpha = 1.0 - min_Gamma - min_Beta;

    // calculate the location of the intersection in world space
    inter.world_intersection = local_to_world_intersection(local_ray, min_t);
    
    // calculate local normal by interpolation in order to calculate 
    // the normal in world space
    Vector3 local_normal = min_v0.normal + min_v1.normal + min_v2.normal;
    inter.world_normal = local_to_world_normal(local_normal);  
    
    inter.local_t = min_t;
    // update the material of the intersection
    inter.material.ambient = material->ambient;
    inter.material.diffuse = material->diffuse;
    inter.material.specular = material->specular;
    inter.material.refractive_index = material->refractive_index;
    
    // calculate the texture color of the intersection 
    Vector2 tex_coord = min_v0.tex_coord * min_Alpha + 
                        min_v1.tex_coord * min_Beta + 
                        min_v2.tex_coord * min_Gamma;
    inter.texture_color = texture_lookup(tex_coord, material);
    
    // update the closest t found 
    *closest_found_t = min_t;
    return; 
}

void Model::render() const
{
    if ( !mesh )
        return;
    if ( material )
        material->set_gl_state();
    mesh->render();
    if ( material )
        material->reset_gl_state();
}

} /* _462 */
