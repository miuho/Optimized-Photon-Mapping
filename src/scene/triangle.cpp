/**
 * @file triangle.cpp
 * @brief Function definitions for the Triangle class.
 *
 * @author Eric Butler (edbutler)
 */

#include "scene/triangle.hpp"
#include "application/opengl.hpp"

namespace _462 {

Triangle::Triangle()
{
    vertices[0].material = 0;
    vertices[1].material = 0;
    vertices[2].material = 0;
}

Triangle::~Triangle() { }

void Triangle::init_bound(Bound& bound) const {
    size_t z;
    // find the minimum and maximum coordinate of the triangle
    for (z = 0; z < 3; z++) {
        Vector3 p = mat.transform_point(vertices[z].position);
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

void Triangle::intersection_test(const Ray& r, Intersection& inter,
                                 real_t *closest_found_t) const {
    // checking the bound first 
    if (intersection_bound_test(r, closest_found_t) == false)
        return;

    // perform the intersection test in object's local space rather than
    // world space, and so the center of the object is <0, 0, 0>
    Ray local_ray = world_to_local_ray(r);
    Vector3 local_e = local_ray.e;
    Vector3 local_d = local_ray.d;
    
    // let the vertices of the triangle be a, b and c, which corresponds
    // to the vertices index as 0, 1, and 2
    // Using Cramer's rule to solve for Beta, Gamma, and t
    real_t a = (vertices[0]).position.x - (vertices[1]).position.x;
    real_t b = (vertices[0]).position.y - (vertices[1]).position.y;
    real_t c = (vertices[0]).position.z - (vertices[1]).position.z;
    real_t d = (vertices[0]).position.x - (vertices[2]).position.x;
    real_t e = (vertices[0]).position.y - (vertices[2]).position.y;
    real_t f = (vertices[0]).position.z - (vertices[2]).position.z;
    real_t g = local_d.x;
    real_t h = local_d.y;
    real_t i = local_d.z;
    real_t j = (vertices[0]).position.x - local_e.x;
    real_t k = (vertices[0]).position.y - local_e.y;
    real_t l = (vertices[0]).position.z - local_e.z;
    real_t M = a * (e * i - h * f) + b * (g * f - d * i) + c * (d * h - e * g);
    
    real_t t = -(f * (a * k - j * b) + e *(j * c - a * l) + 
               d * (b * l - k * c))/M;
    // check if we already found a closer t
    if (t < t0 || t > *closest_found_t) {
        return;
    } 

    real_t Gamma = (i * (a * k - j * b) + h * (j * c - a * l) + 
                   g * (b * l - k * c))/M;
    if (Gamma < 0.0 || Gamma > 1.0) {
        return;
    }

    real_t Beta = (j * (e * i - h * f) + k * (g * f - d * i) + 
                  l * (d * h - e * g))/M;
    if (Beta < 0.0 || Beta > (1.0 - Gamma)) {
        return;
    }
    real_t Alpha = 1.0 - Beta - Gamma;
     
    // calculate the location of the intersection in world space
    inter.world_intersection = local_to_world_intersection(local_ray, t);
    
    // calculate local normal by interpolation in order to calculate 
    // the normal in world space
    Vector3 local_normal = (vertices[0]).normal + 
                           (vertices[1]).normal + 
                           (vertices[2]).normal;
    inter.world_normal = local_to_world_normal(local_normal);
    
    inter.local_t = t;
    // update the material of the intersection by interpolation
    inter.material.ambient = 
        (vertices[0]).material->ambient * Alpha +
        (vertices[1]).material->ambient * Beta +
        (vertices[2]).material->ambient * Gamma;

    inter.material.diffuse = 
        (vertices[0]).material->diffuse * Alpha +
        (vertices[1]).material->diffuse * Beta +
        (vertices[2]).material->diffuse * Gamma;

    inter.material.specular =
        (vertices[0]).material->specular * Alpha +
        (vertices[1]).material->specular * Beta +
        (vertices[2]).material->specular * Gamma;


    inter.material.refractive_index = 
        (vertices[0]).material->refractive_index * Alpha +
        (vertices[1]).material->refractive_index * Beta +
        (vertices[2]).material->refractive_index * Gamma;

    
    // calculate the texture color of the intersection 
    Vector2 tex_coord = (vertices[0]).tex_coord * Alpha + 
                        (vertices[1]).tex_coord * Beta + 
                        (vertices[2]).tex_coord * Gamma;
    inter.texture_color = Alpha * texture_lookup(tex_coord, 
                                                 (vertices[0]).material);
    inter.texture_color += Beta * texture_lookup(tex_coord, 
                                                 (vertices[1]).material);
    inter.texture_color += Gamma * texture_lookup(tex_coord, 
                                                  (vertices[2]).material);

    // update the closest t found
    *closest_found_t = t;
    return;
}

void Triangle::render() const
{
    bool materials_nonnull = true;
    for ( int i = 0; i < 3; ++i )
        materials_nonnull = materials_nonnull && vertices[i].material;

    // this doesn't interpolate materials. Ah well.
    if ( materials_nonnull )
        vertices[0].material->set_gl_state();

    glBegin(GL_TRIANGLES);

    glNormal3dv( &vertices[0].normal.x );
    glTexCoord2dv( &vertices[0].tex_coord.x );
    glVertex3dv( &vertices[0].position.x );

    glNormal3dv( &vertices[1].normal.x );
    glTexCoord2dv( &vertices[1].tex_coord.x );
    glVertex3dv( &vertices[1].position.x);

    glNormal3dv( &vertices[2].normal.x );
    glTexCoord2dv( &vertices[2].tex_coord.x );
    glVertex3dv( &vertices[2].position.x);

    glEnd();

    if ( materials_nonnull )
        vertices[0].material->reset_gl_state();
}

} /* _462 */
