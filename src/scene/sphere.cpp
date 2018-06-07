/**
 * @file sphere.cpp
 * @brief Function defnitions for the Sphere class.
 *
 * @author Kristin Siu (kasiu)
 * @author Eric Butler (edbutler)
 */

#include "scene/sphere.hpp"
#include "application/opengl.hpp"

namespace _462 {

#define SPHERE_NUM_LAT 80
#define SPHERE_NUM_LON 100

#define SPHERE_NUM_VERTICES ( ( SPHERE_NUM_LAT + 1 ) * ( SPHERE_NUM_LON + 1 ) )
#define SPHERE_NUM_INDICES ( 6 * SPHERE_NUM_LAT * SPHERE_NUM_LON )
// index of the x,y sphere where x is lat and y is lon
#define SINDEX(x,y) ((x) * (SPHERE_NUM_LON + 1) + (y))
#define VERTEX_SIZE 8
#define TCOORD_OFFSET 0
#define NORMAL_OFFSET 2
#define VERTEX_OFFSET 5

static unsigned int Indices[SPHERE_NUM_INDICES];
static float Vertices[VERTEX_SIZE * SPHERE_NUM_VERTICES];

static void init_sphere()
{
    static bool initialized = false;
    if ( initialized )
        return;

    for ( int i = 0; i <= SPHERE_NUM_LAT; i++ ) {
        for ( int j = 0; j <= SPHERE_NUM_LON; j++ ) {
            real_t lat = real_t( i ) / SPHERE_NUM_LAT;
            real_t lon = real_t( j ) / SPHERE_NUM_LON;
            float* vptr = &Vertices[VERTEX_SIZE * SINDEX(i,j)];

            vptr[TCOORD_OFFSET + 0] = lon;
            vptr[TCOORD_OFFSET + 1] = 1-lat;

            lat *= PI;
            lon *= 2 * PI;
            real_t sinlat = sin( lat );

            vptr[NORMAL_OFFSET + 0] = vptr[VERTEX_OFFSET + 0] = sinlat * sin( lon );
            vptr[NORMAL_OFFSET + 1] = vptr[VERTEX_OFFSET + 1] = cos( lat ),
            vptr[NORMAL_OFFSET + 2] = vptr[VERTEX_OFFSET + 2] = sinlat * cos( lon );
        }
    }

    for ( int i = 0; i < SPHERE_NUM_LAT; i++ ) {
        for ( int j = 0; j < SPHERE_NUM_LON; j++ ) {
            unsigned int* iptr = &Indices[6 * ( SPHERE_NUM_LON * i + j )];

            unsigned int i00 = SINDEX(i,  j  );
            unsigned int i10 = SINDEX(i+1,j  );
            unsigned int i11 = SINDEX(i+1,j+1);
            unsigned int i01 = SINDEX(i,  j+1);

            iptr[0] = i00;
            iptr[1] = i10;
            iptr[2] = i11;
            iptr[3] = i11;
            iptr[4] = i01;
            iptr[5] = i00;
        }
    }

    initialized = true;
}

Sphere::Sphere()
    : radius(0), material(0) {}

Sphere::~Sphere() {}

void Sphere::init_bound(Bound& bound) const {
    size_t z;
    real_t a, b, c;
    // the cube bounding box that encloses the sphere has 8 vertices which
    // are all possible linear combinations of <radius, radius, radius>
    for (z = 0; z < 8; z++) {
        a = (z < 4) ? 1.0 : -1.0; 
        b = (z % 4 <= 1) ? 1.0 : -1.0; 
        c = (z % 2 == 0) ? 1.0 : -1.0; 
        Vector3 q = Vector3(a * radius, b * radius, c * radius);
        Vector3 p = mat.transform_point(q);
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

void Sphere::intersection_test(const Ray& r, Intersection& inter,
                               real_t *closest_found_t) const {
    // checking the bound first 
    if (intersection_bound_test(r, closest_found_t) == false)
        return;
    
    real_t t;
    // perform the intersection test in object's local space rather than
    // world space, and so the center of the object is <0, 0, 0>
    Ray local_ray = world_to_local_ray(r);
    Vector3 local_e = local_ray.e;
    Vector3 local_d = local_ray.d;

    real_t A = dot(local_d, local_d);
    real_t B = 2.0 * dot(local_d, local_e);
    real_t C = dot(local_e, local_e) - radius * radius;
    // calculate the discriminant to check for intersection
    real_t discriminant = B * B - 4.0 * A * C;

    // sphere do not intersect if discriminant is negative
    if (discriminant < 0.0) {
        return;
    }
    // ray grazes the sphere, exactly 1 point if discrimiant is zero
    else if (discriminant == 0.0) {
        t = (-1.0 * B) / (2.0 * A);
        
        // return if the t is further than what we found already
        if (t < t0 || t > *closest_found_t) {
            return;
        }
        
        // calculate the location of the intersectino in world space
        inter.world_intersection = local_to_world_intersection(local_ray, t);
    }
    // 2 intersections (one entering and one leaving the sphere) if 
    // discrimiant is positive
    else {
        real_t t1 = ((-1.0 * B) + sqrt(discriminant)) / (2.0 * A);
        real_t t2 = ((-1.0 * B) - sqrt(discriminant)) / (2.0 * A);
        
        if (t1 < t0 || t1 > *closest_found_t) {
            // t1 is out of bound
            if (t2 < t0 || t2 > *closest_found_t) {
                // t1 and t2 are out of bound
                return;
            }
            else {
                // only t2 is in bound
                t = t2;
            }
        }       
        else {
            // t1 is in bound
            if (t2 < t0 || t2 > *closest_found_t) {
                // t2 is out of bound
                t = t1;
            }
            else {
                // t1 and t2 are in bound
                t = (t1 < t2) ? t1 : t2;
            }
        } 
        
        // calculate the location of the intersection in world space
        inter.world_intersection = local_to_world_intersection(local_ray, t);
    }
        
    // local normal is the direction pointing away from the center of sphere
    // which is the local intersection minus <0, 0, 0>, and so we use the
    // local normal to calculte the normal in world space
    Vector3 local_normal = local_e + t * local_d;
    inter.world_normal = local_to_world_normal(local_normal);
    
    inter.local_t = t;
    // update the material of the intersection
    inter.material.ambient = material->ambient;
    inter.material.diffuse = material->diffuse;
    inter.material.specular = material->specular;
    inter.material.refractive_index = material->refractive_index;
        
    // calculate the texture color of the sphere
    Vector2 tex_coord = Vector2();
    real_t Theta = acos(local_normal.z / radius);
    real_t Phi = atan2(local_normal.y, local_normal.x);
    // add 2pi to Phi if it is negative
    if (Phi < 0) 
        Phi += 2.0 * PI;
    tex_coord.x = Phi/(2.0 * PI);
    tex_coord.y = (PI - Theta)/ PI;
    inter.texture_color = texture_lookup(tex_coord, material); 

    // update the closest t found
    *closest_found_t = t;
    return;
}

void Sphere::render() const
{
    // create geometry if we haven't already
    init_sphere();

    if ( material )
        material->set_gl_state();

    // just scale by radius and draw unit sphere
    glPushMatrix();
    glScaled( radius, radius, radius );
    glInterleavedArrays( GL_T2F_N3F_V3F, VERTEX_SIZE * sizeof Vertices[0], Vertices );
    glDrawElements( GL_TRIANGLES, SPHERE_NUM_INDICES, GL_UNSIGNED_INT, Indices );
    glPopMatrix();

    if ( material )
        material->reset_gl_state();
}

} /* _462 */

