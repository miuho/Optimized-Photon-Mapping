/**
 * @file raytacer.cpp
 * @brief Raytracer class
 *
 * 
 *
 * @author HingOn Miu (hmiu)
 * 
 */

#include "raytracer.hpp"
#include "scene/scene.hpp"

#include <SDL_timer.h>
#include <iostream>
#include <random>

#ifdef OPENMP // just a defense in case OpenMP is not installed.

#include <omp.h>

#endif
namespace _462 {

// max number of threads OpenMP can use. Change this if you like.
#define MAX_THREADS 16
// the number of shadow ray to sample the spherical light
#define NUM_SHADOWRAY 5
// split the photon array according to x, y, or z axis
#define X_AXIS 0
#define Y_AXIS 1
#define Z_AXIS 2
// indicate the leaf node
#define LEAF 3
// set the number of photons fire from light source
#define NUM_PHOTON 100000
// set the scaling factor for photon internsity
#define FACTOR 0.001
// set a upper bound for square distance calculation
#define BIG_NUM 999999.9
// the disk range is about a 2 degree offset from 90 degree (88 - 92)
#define EPSILON 0.033
// set the number of photon neighbors to use for radiance estimate
#define NUM_PHOTON_RADIANCE 500
// differentiate the previous bounce for photon tracing
#define DIFFUSE 15462
#define SPECULAR 15662

static const unsigned STEP_SIZE = 8;
// map to store caustic photons
std::vector<Photon> caustic_map;
// map to store global illumination photons
std::vector<Photon> global_illum_map;
// records whether the photon maps are initialized
static bool maps_initialized = false;
static void emit_photons(const Scene* scene);

Raytracer::Raytracer()
    : scene(0), width(0), height(0) { }

// random real_t in [0, 1)
static inline real_t random()
{
    return real_t(rand())/RAND_MAX;
}

/**
 * Find the maximum of 2 inputs.
 * @param a The first input.
 * @param b The second input.
 * @return maximum of 3 inputs
 */
inline real_t max2(real_t a, real_t b) {
    return (a > b) ? a : b;
}

/**
 * Find the maximum of 3 inputs.
 * @param a The first input.
 * @param b The second input.
 * @param c The third input.
 * @return maximum of 3 inputs
 */
inline real_t max3(real_t a, real_t b, real_t c) {
    return max2(max2(a, b), c);
}

Raytracer::~Raytracer() { }

/**
 * Initializes the raytracer for the given scene. Overrides any previous
 * initializations. May be invoked before a previous raytrace completes.
 * @param scene The scene to raytrace.
 * @param width The width of the image being raytraced.
 * @param height The height of the image being raytraced.
 * @return true on success, false on error. The raytrace will abort if
 *  false is returned.
 */
bool Raytracer::initialize(Scene* scene, size_t num_samples,
			   size_t width, size_t height)
{
    /*
     * omp_set_num_threads sets the maximum number of threads OpenMP will
     * use at once.
     */
#ifdef OPENMP
    omp_set_num_threads(MAX_THREADS);
#endif
    this->scene = scene;
    this->num_samples = num_samples;
    this->width = width;
    this->height = height;

    current_row = 0;

    Ray::init(scene->camera);
    scene->initialize();

    // TODO any initialization or precompuation before the trace
    if (maps_initialized == false) {
	// initializes the photon maps
	caustic_map.reserve(NUM_PHOTON * MAX_DEPTH);
	global_illum_map.reserve(NUM_PHOTON * MAX_DEPTH);
	 
	// do photon tracing for light source
	emit_photons(scene);
	maps_initialized = true;
    }
    
    return true;
}

/**
 * Compare the x coordinates of two photons.
 * @param i The first photon.
 * @param j The second photon.
 * @return true if i smaller than j, false otherwise.
 */
inline bool compare_x(Photon i, Photon j) {
    return i.p.x < j.p.x; 
}

/**
 * Compare the y coordinates of two photons.
 * @param i The first photon.
 * @param j The second photon.
 * @return true if i smaller than j, false otherwise.
 */
inline bool compare_y(Photon i, Photon j) {
    return i.p.y < j.p.y; 
}

/**
 * Compare the z coordinates of two photons.
 * @param i The first photon.
 * @param j The second photon.
 * @return true if i smaller than j, false otherwise.
 */
inline bool compare_z(Photon i, Photon j) {
    return i.p.z < j.p.z; 
}

/**
 * Balance the kd tree.
 * @param L The photon map.
 * @param begin The begin index of the array.
 * @param end The end index of the array.
 * @return void.
 */
void construct_kdtree(std::vector<Photon>& L, unsigned begin, unsigned end) {
    if (end - begin == 0) {
        return;
    }
    if (end - begin == 1) {
        // indicate the leaf node
        L[begin].flag = LEAF;
        return;
    }
    
    // calculate the variance
    unsigned median = begin + (end - begin)/2;
    real_t x_avg = 0.0, y_avg = 0.0, z_avg = 0.0;
    real_t x_var = 0.0, y_var = 0.0, z_var = 0.0;
    real_t n = end - begin;
    std::vector<Photon>::iterator a = L.begin() + begin;
    std::vector<Photon>::iterator b = L.begin() + end;
    std::vector<Photon>::iterator it;
    for (it = a; it != b; ++it) {
        x_avg += (*it).p.x;
        y_avg += (*it).p.y;
        z_avg += (*it).p.z;
    }
    for (it = a; it != b; ++it) {
        x_var += ((*it).p.x - x_avg) * ((*it).p.x - x_avg);
        y_var += ((*it).p.y - y_avg) * ((*it).p.y - y_avg);
        z_var += ((*it).p.z - z_avg) * ((*it).p.z - z_avg);
    }
    x_var /= n; 
    y_var /= n; 
    z_var /= n; 
    
    // find the dimension with maximum variance 
    real_t max_var = max3(x_var, y_var, z_var);

    // split the dimension and indicate the splitting axis
    if (max_var == x_var) {
        std::sort(L.begin() + begin, L.begin() + end, compare_x);
        L[median].flag = X_AXIS;
    }
    if (max_var == y_var) {
        std::sort(L.begin() + begin, L.begin() + end, compare_y);
        L[median].flag = Y_AXIS;
    }
    if (max_var == z_var) {
        std::sort(L.begin() + begin, L.begin() + end, compare_z);
        L[median].flag = Z_AXIS;
    }

    // recurse on left and right children 
    construct_kdtree(L, begin, median);        
    construct_kdtree(L, median + 1, end);  
    return; 
}

/**
 * Swap two elements in a max heap.
 * @param neighbors The max heap array.
 * @param a The index of first element.
 * @param b The index of second element.
 * @return void.
 */
inline void heap_swap(Neighbor *neighbors, int a, int b) {
    unsigned a_i = (neighbors[a]).i;
    real_t a_s = (neighbors[a]).sq_dis;
    
    (neighbors[a]).i = (neighbors[b]).i;
    (neighbors[a]).sq_dis = (neighbors[b]).sq_dis;
    (neighbors[b]).i = a_i;
    (neighbors[b]).sq_dis = a_s;  
}

/**
 * Remove an element in a max heap.
 * @param neighbors The max heap array.
 * @param size The size of the max heap.
 * @return void.
 */
void heap_remove(Neighbor *neighbors, int *size) {
    // move the last element to the root node so that
    // the max element is replaced
    (neighbors[0]).i = (neighbors[*size - 1]).i;
    (neighbors[0]).sq_dis = (neighbors[*size - 1]).sq_dis;
    *size = *size - 1; 
    
    int i = 0;
    int left, right, bigger;
    real_t i_val, left_val, right_val;
    // swap the root node element downward until it has no children
    // or both children have smaller values
    while (1) {
        left = 2*i + 1;
        right = 2*i + 2;
        if (left >= *size && right >= *size) {
            // i is a leaf node (has no child)
            return;
        }
        i_val = (neighbors[i]).sq_dis;
        left_val = (left < *size) ? (neighbors[left]).sq_dis : -1.0;
        right_val = (right < *size) ? (neighbors[right]).sq_dis : -1.0; 
        if (i_val >= left_val && i_val >= right_val) {
            // i is bigger than both children
            return;
        } 
        if (left_val == -1.0 && right_val != -1.0) {
            // i is smaller than right child 
            heap_swap(neighbors, i, right);
            i = right;
        } 
        if (left_val != -1.0 && right_val == -1.0) {
            // i is smaller than left child 
            heap_swap(neighbors, i, left);
            i = left;
        }
        else {
            // i is smaller than at least one of the child
            bigger = (left_val > right_val) ? left : right;
            heap_swap(neighbors, i, bigger);
            i = bigger; 
        }
    }
         
    return;
}

/**
 * Insert an element in a max heap.
 * @param neighbors The max heap array.
 * @param size The size of the max heap.
 * @param e The index of photon in map.
 * @param e_dis The square distance of photon and intersection.
 * @return void.
 */
inline void heap_add(Neighbor *neighbors, int *size,
              unsigned e, real_t e_dis) {
    // insert a new element to the last element of max heap 
    int i = *size;
    (neighbors[i]).i = e;
    (neighbors[i]).sq_dis = e_dis;
    *size = *size + 1;
    
    int parent;
    real_t i_val, parent_val;
    // swap the last element upward unitl it reaches the root node
    // or has a larger parent
    while (1) {
        if (i == 0) {
            // reached root node
            return;
        }
        parent = (i - 1)/2;
        i_val = (neighbors[i]).sq_dis;
        parent_val = (neighbors[parent]).sq_dis;
        if (parent_val >= i_val) {
            // parent is bigger than i
            return;
        }
        
        heap_swap(neighbors, i, parent);
        i = parent;
    }

    return;
}

/**
 * Add a new neighbor to the nearest neighbor max heap.
 * @param p The interection position.
 * @param norm The interection normal.
 * @param neighbors The max heap array.
 * @param L The photon map.
 * @param e The photon index of the map.
 * @param D The maximum square distance of photon in neighor heap.
 * @param size The size of the max heap.
 * @return void.
 */
inline void add_neighbor(Vector3 p, Vector3 norm, Neighbor *neighbors, 
                  std::vector<Photon>& L, unsigned e, real_t *D, int *size) {
    // disk check
    if (abs(dot(norm, normalize(p -  L[e].p))) > EPSILON) {
        return;
    }
    
    real_t e_dis = squared_distance(p, (L[e]).p);
    if (*size < NUM_PHOTON_RADIANCE  || e_dis < *D) {
        // maintain the size of the max heap 
        if (*size == NUM_PHOTON_RADIANCE) {
            heap_remove(neighbors, size); 
        }
    
        heap_add(neighbors, size, e, e_dis);
        
        // update the maximum square distance
        *D = (neighbors[0]).sq_dis;
    }
}

/**
 * Get the split value.
 * @param L The photon map.
 * @param i The photon index of the map.
 * @param axis The splitting axis.
 * @return The split value.
 */
inline real_t get_split(std::vector<Photon>& L, unsigned i, int axis) {
    if (axis == X_AXIS)
        return (L[i]).p.x;
    if (axis == Y_AXIS)
        return (L[i]).p.y;
    if (axis == Z_AXIS)
        return (L[i]).p.z;
    return 0.0;
}

/**
 * Get the coordinate of intersection of the splitting dimension.
 * @param p The interection position.
 * @param axis The splitting axis.
 * @return The coordinate of the splitting dimension.
 */
inline real_t get_p(Vector3 p, int axis) {
    if (axis == X_AXIS)
        return p.x;
    if (axis == Y_AXIS)
        return p.y;
    if (axis == Z_AXIS)
        return p.z;
    return 0.0;
}

/**
 * Get the nearest photon neighbors of a position.
 * @param p The interection position.
 * @param norm The interection normal.
 * @param neighbors The max heap array.
 * @param L The photon map.
 * @param e The photon index of the map.
 * @param D The maximum square distance of photon in neighor heap.
 * @param size The size of the max heap.
 * @return void.
 */
void lookup(Vector3 p, Vector3 norm, Neighbor *neighbors, std::vector<Photon>& L, 
            unsigned begin, unsigned end, real_t *D, int *size) {
    if (begin == end)
        return;
    if (begin + 1 == end)
	// add photon at leaf node to neighbors heap
        return add_neighbor(p, norm, neighbors, L, begin, D, size);
    
    unsigned median = begin + (end - begin)/2;
    // get splitting axis
    int flag = (L[median]).flag;
    real_t split_value = get_split(L, median, flag);
    real_t p_value = get_p(p, flag);
    // check which side of the splitting axis to traverse
    if (p_value <= split_value) {
	// traverse left sub-tree first
        lookup(p, norm, neighbors, L, begin, median, D, size);
	// add the current node
        add_neighbor(p, norm, neighbors, L, median, D, size);   
        // return if neighbors heap is full and all nodes in the
	// right sub-tree is further than those in neighbors heap
        if (*size >= NUM_PHOTON_RADIANCE && 
            (p_value - split_value) * (p_value - split_value) > *D) {
            return;
        }
        // traverse right sub-tree
        return lookup(p, norm, neighbors, L, median + 1, end, D, size);
    }
    else {
	// traverse right sub-tree first
        lookup(p, norm, neighbors, L, median + 1, end, D, size);
	// add the current node
        add_neighbor(p, norm, neighbors, L, median, D, size);
        // return if neighbors heap is full and all nodes in the
	// left sub-tree is further than those in neighbors heap
        if (*size >= NUM_PHOTON_RADIANCE && 
            (p_value - split_value) * (p_value - split_value) > *D) {
            return;
        }
        // traverse left sub-tree
        return lookup(p, norm, neighbors, L, begin, median, D, size);   
    }
}

/**
 * Get the color from global illumination.
 * @param p The interection position.
 * @param norm The interection normal.
 * @param neighbors The max heap array.
 * @return The global illumination color at intersection.
 */
Color3 global_illumination(Neighbor *neighbors, Vector3 p, Vector3 norm) {
    Color3 result = Color3::Black();
    // D is the furthest photon squared distance to p in neighbors
    real_t D = BIG_NUM;
    int size = 0;
    // find the nearest photons
    lookup(p, norm, neighbors, global_illum_map, 
           0, global_illum_map.size(), &D, &size);
    if (size == 0) 
        return result;
    for (int i = 0; i < size; i++) {
        Photon ph = global_illum_map[(neighbors[i]).i];
        result += ph.c * max2(dot(norm, ph.dir), 0.0);
    }
    return result * (1.0 / (D * PI));
}

/**
 * Get the color from caustic.
 * @param p The interection position.
 * @param norm The interection normal.
 * @param neighbors The max heap array.
 * @return The caustic color at intersection.
 */
Color3 caustic(Neighbor *neighbors, Vector3 p, Vector3 norm) {
    Color3 result = Color3::Black();
    // D is the furthest photon squared distance to p in neighbors
    real_t D = BIG_NUM;
    int size = 0;
    // find the nearest photons
    lookup(p, norm, neighbors, caustic_map, 
           0, caustic_map.size(), &D, &size);
    if (size == 0) 
        return result;
    for (int i = 0; i < size; i++) {
        Photon ph = caustic_map[(neighbors[i]).i];
        result += ph.c * max2(dot(norm, ph.dir), 0.0);
    }
    return result * (1.0 / (D * PI));
}

/**
 * Cast a ray to the scene and check if any object is hit.
 * @param scene The scene to raytrace.
 * @param r The ray to cast.
 * @param inter The details of intersected object.
 * @return true if hit an object, and false otherwise.
 */
inline bool cast_ray(const Scene* scene, const Ray& r, Intersection& inter) {
    size_t i;
    // record the closest t found
    real_t closest_found_t = t1;
    // check if the ray does hit an object 
    for (i = 0; i < scene->num_geometries(); i++) {
        // instead of comparing distance here inefficiently,
        // intersection_test check closes_found_t before editing
        // inter, and so inter will be the closest intersection
        // after the for loop.
        ((scene->get_geometries())[i])->
            intersection_test(r, inter, &closest_found_t);
    }
    
    // check if the closest_found_t is still the upper time bound
    // to confirm whether there is an intersection
    return closest_found_t > t0 && closest_found_t < t1;
}

/**
 * Use Monte Carlo method to sample the spherical light.
 * @param light The spherical light.
 * @return The new light position.
 */
inline Vector3 monte_carlo_sampling(const SphereLight light) {
    real_t x = random_gaussian();
    real_t y = random_gaussian();
    real_t z = random_gaussian();
    // generate a random point on a unit sphere
    Vector3 surface_p = normalize(Vector3(x,y,z));
    // scale the point by sphere radius and transform it by sphere center 
    return surface_p * light.radius + light.position;   
}

/**
 * Use Blinn-Phong illumination to calculate direct illumination.
 * @param scene The scene to render.
 * @param r The ray that hits an object.
 * @param inter The details of the intersected object.
 * @return The color of the intersection.
 */
Color3 direct_illumination(const Scene *scene, const Ray &r, 
                           Intersection& inter) {
    // support the basic Blinn-Phong illumination
    // compute the ambient term
    Color3 a = scene->ambient_light * inter.material.ambient;
    
    size_t i, j;
    real_t dis, non_blocked;
    Vector3 dir, volume_light_dir;
    Color3 d = Color3::Black();
    Color3 light_color;
    real_t N_L;
    Intersection tmp;
    Ray light_r;
    // compute the diffuse term
    for (i = 0; i < scene->num_lights(); i++) {
        // light direction
        dir = scene->get_lights()[i].position - 
              inter.world_intersection;
        // distance from the light
        dis = length(dir);
	
	// record the number of non-blocked shadow rays
        non_blocked = NUM_SHADOWRAY;
        for (j = 0; j < NUM_SHADOWRAY; j++) {     
            // use monte carlo method to sample volumne light source position
            volume_light_dir = monte_carlo_sampling(scene->get_lights()[i])
                               - inter.world_intersection; 
            light_r = Ray(inter.world_intersection, volume_light_dir); 
            // check the shadow, whether the light is blocked by closer object
            if (cast_ray(scene, light_r, tmp) == true &&
                dis > length(tmp.world_intersection-inter.world_intersection)) {
                non_blocked--;
            }
        } 
        
        // compute the color of the light with attenuation
        light_color = (scene->get_lights()[i]).color *  
             (1.0 / (scene->get_lights()[i].attenuation.constant +
             scene->get_lights()[i].attenuation.linear * dis +
             scene->get_lights()[i].attenuation.quadratic * dis * dis));       
        // dot product of intersection normal and light vector
        N_L = dot(inter.world_normal, normalize(dir));
                  //inter.world_intersection - 
        // max(N_L, 0.0)
        if (N_L < 0.0) {
            N_L = 0.0;
        }        
        
        d += (non_blocked/(real_t)NUM_SHADOWRAY) * light_color * 
             inter.material.diffuse * N_L;
    }
    
    return inter.texture_color * (a + d);
}

/**
 * Calculate the reflected ray.
 * @param e The intersection location.
 * @param d The direction of incoming ray.
 * @param n The normal vector on surface.
 * @return The refected ray.
 */
inline Ray reflect(Vector3 e, Vector3 d, Vector3 n) {
    // reflected_dir is d - 2*(d.n)*n
    Vector3 reflect_dir = d - 2.0 * dot(d, n) * n;
    return Ray(e, reflect_dir);
}

/**
 * Check whether there is refracted ray from the intersected object.
 * @param d The direction of incoming ray.
 * @param norm The normal vector of surface.
 * @param n The current medium's refractive index.
 * @param n The intersected medium's refractive index.
 * @param t The refracted ray direction.
 * @return true if there is refracted ray, and false otherwise.
 */
inline bool refract(Vector3 d, Vector3 norm, real_t n, real_t nt, Vector3* t) {
    // if the number under the square root is negative, there is no
    // refracted ray and all of the energy is reflected, which is
    // known as total internal reflection.
    real_t d_n = dot(d, norm);
    real_t sq_rt = 1.0 - (n * n * (1.0 - 
                   (d_n * d_n)))/ (nt * nt);

    if (sq_rt < 0.0) {
        return false;
    }
    else {
        Vector3 tmp = 
            (n * (d - norm * d_n)) / nt - norm * sqrt(sq_rt);
        t->x = tmp.x;
        t->y = tmp.y;
        t->z = tmp.z;
        return true;
    }
}

/**
 * Get the pixel color by calculating direct illumination, reflection,
 * and refraction.
 * @param scene The scene to render.
 * @param r The ray to cast.
 * @param depth The recursion depth.
 * @return The color of the pixel.
 */
Color3 get_pixel_color(const Scene *scene, const Ray& r, int depth, Neighbor *neighbors) {
    // the maximum recursive depth for tracing the ray from screen
    if (depth > MAX_DEPTH) {
         return Color3::Black();
    }
    
    Intersection inter; 
    // nothing is hit by the eye ray 
    if (cast_ray(scene, r, inter) == false) {
        return scene->background_color;
    }
    
    // compute the reflected ray 
    Ray reflect_r = reflect(inter.world_intersection, r.d, 
                            inter.world_normal); 

    // handle opaque objects
    if (inter.material.refractive_index == 0.0) {
        // compute the color by direct illumination
        //Color3 direct = Color3::Black();
        Color3 direct = direct_illumination(scene, r, inter);
        // compute the color by global illumination
        //Color3 global = Color3::Black();
        Color3 global = global_illumination(neighbors, 
            inter.world_intersection, inter.world_normal);
        // compute the color by caustics
        //Color3 caus = Color3::Black();
        Color3 caus = caustic(neighbors, 
            inter.world_intersection, inter.world_normal);
        // sum the color from direct illumination and 
        // color from specular term's reflected ray computation
        return direct + inter.texture_color * inter.material.diffuse * 
            (global + caus) + inter.texture_color * inter.material.specular *
            get_pixel_color(scene, reflect_r, depth + 1, neighbors);
    }
    // handle dielectrics
    else {
        real_t c, n, nt;
        Vector3 t;
        bool entering;
        
        // entering the dielectric 
        if (dot(r.d, inter.world_normal) < 0.0) {
            entering = true;
            n = scene->refractive_index;
            nt = inter.material.refractive_index;
            refract(r.d, inter.world_normal, n, nt, &t);
            c = dot(-r.d, inter.world_normal);

        }
        // leaving the dialectric
        else {
            entering = false;
            n = inter.material.refractive_index;
            nt = scene->refractive_index;
            if (refract(r.d, -inter.world_normal, n, nt, &t)) {
                c = dot(t, inter.world_normal); 
            }
            // total internal reflection
            else {
                // R = 1 for total internal reflection
                return inter.texture_color * inter.material.specular *
                       get_pixel_color(scene, reflect_r, depth + 1, neighbors);
            }
        }
        
        // compute Fresnel coefficient R
        real_t R0 = (real_t)pow((double)((nt - 1.0)/(nt + 1.0)), 2.0);
        real_t R = R0 + (1.0 - R0) * (real_t)pow(1.0 - c, 5.0);
        Ray refract_r = Ray(inter.world_intersection, t);
        
        // schlick approximation of the Fresnel effect
        // point color = R * reflection color + (1-R) * refraction color
        if (entering) {
            // entering object means the new current refractive index
            // would be the one next on stack
            return inter.texture_color * inter.material.specular * ( 
                R * get_pixel_color(scene, reflect_r, depth + 1, neighbors) +
                (1.0 - R) * get_pixel_color(scene, refract_r, depth + 1, neighbors));
        }
        else {
            // leaving object means the new current refractive index
            // would be the one previous on stack
            return inter.texture_color * inter.material.specular * (
                R * get_pixel_color(scene, reflect_r, depth + 1, neighbors) +
                (1.0 - R) * get_pixel_color(scene, refract_r, depth + 1, neighbors));
        }

    }
}

/**
 * Generate a photon ray.
 * @param light The spherical light in scene.
 * @return The photon ray.
 */
inline Ray photon_ray(const SphereLight light) {
    Vector3 p;
    // a random position on unit sphere
    p.x = 2.0 * random() - 1.0;
    p.y = 2.0 * random() - 1.0;
    p.z = 2.0 * random() - 1.0;
    p = normalize(p);
    // initial ray for a photon
    return Ray(p + light.position, p);
}

/**
 * Random unit direction on a hemisphere.
 * @param normal The hemisphere direction.
 * @return The unit direction.
 */
inline Vector3 uniformSampleHemiSphere(const Vector3& normal) {
    Vector3 p;
    p.x = 2.0 * random() - 1.0;
    p.y = 2.0 * random() - 1.0;
    p.z = 2.0 * random() - 1.0;
    p = normalize(p);   
    if (dot(p, normal) < 0.0)
        p = -p;
    return p;
}

/**
 * Calculate the probability of reflection.
 * @param coef Either the diffuse or specular material color.
 * @param power The intensity of the photon.
 * @return The reflection probability.
 */
inline real_t get_reflect_probability(Color3 coef, Color3 power) {
    real_t top = max3(coef.r * power.r, coef.g * power.g, coef.b * power.b);
    real_t bot = max3(power.r, power.g, power.b);
    return top/bot;
}

/**
 * Trace the photon ray and insert the photons into photon maps
 * at diffuse surfaces.
 * @param scene The scene to render.
 * @param r The ray to cast.
 * @param ph The photon.
 * @param depth The recursion depth.
 * @param previous_bounce The previous bounce of photon.
 * @return void.
 */
void trace_photon(const Scene *scene, const Ray& r, 
                  Photon& ph, int depth, int previous_bounce) {
    // the maximum recursive depth for tracing the ray from screen
    if (depth > MAX_DEPTH) {
         return;
    }
    
    Intersection inter; 
    // nothing is hit by the eye ray 
    if (cast_ray(scene, r, inter) == false) {
        return;
    }
    // find the photon's position
    ph.p = inter.world_intersection;

    // handle opaque objects
    if (inter.material.refractive_index == 0.0) {
        real_t P_d = get_reflect_probability(inter.material.diffuse,
                                             ph.c);
        real_t P_s = get_reflect_probability(inter.material.specular,
                                             ph.c);
        real_t x = random();
        Ray reflect_r;
        Photon reflected_ph;
        int bounce;
        bool absorb = false;
        if (x >= 0 && x < P_d) {
            // [0, P_d] diffuse reflection
	    // reflected photon ray has random direction
            reflect_r = Ray(inter.world_intersection, 
                            uniformSampleHemiSphere(inter.world_normal));
            reflected_ph = Photon(ph.c * inter.material.diffuse *
                                  inter.texture_color);
            bounce = DIFFUSE; 
        }   
        if (x >= P_d && x < (P_s + P_d)) {
            // [P_d, P_s + P_d] specular reflection
            reflect_r = reflect(inter.world_intersection, normalize(r.d), 
                                inter.world_normal); 
            reflected_ph = Photon(ph.c * inter.material.specular * 
                                  inter.texture_color);
            bounce = SPECULAR;   
        }
        if (x >= (P_s + P_d)) { 
            // [P_s + P_d, 1] absorption
            absorb = true;
        }
        
        // do not add to map for first intersection 
        if (depth != 0) {
            ph.dir = -normalize(r.d);
            if (previous_bounce == DIFFUSE)
                global_illum_map.push_back(ph);
            if (previous_bounce == SPECULAR)
                caustic_map.push_back(ph);
        }
        
	//terminate tracing if the photon is absorbed
        if (absorb)
            return;
        else  
            return trace_photon(scene, reflect_r, reflected_ph, depth + 1, bounce);
    }
    // handle dielectrics
    else {
        real_t c, n, nt;
        Vector3 t;
        bool entering;
         
        // compute the reflected ray 
        Ray reflect_r = reflect(inter.world_intersection, normalize(r.d), 
                                inter.world_normal);     
        // adjust photon power
        ph.c = ph.c * inter.texture_color * inter.material.specular;

        // entering the dielectric 
        if (dot(r.d, inter.world_normal) < 0.0) {
            entering = true;
            n = scene->refractive_index;
            nt = inter.material.refractive_index;
            refract(normalize(r.d), inter.world_normal, n, nt, &t);
            c = dot(-normalize(r.d), inter.world_normal);

        }
        // leaving the dialectric
        else {
            entering = false;
            n = inter.material.refractive_index;
            nt = scene->refractive_index;
            if (refract(normalize(r.d), -inter.world_normal, n, nt, &t)) {
                c = dot(normalize(t), inter.world_normal); 
            }
            // total internal reflection
            else {
                // R = 1 for total internal reflection
                return trace_photon(scene, reflect_r, ph, depth + 1, SPECULAR);
            }
        }
        
        // compute Fresnel coefficient R
        real_t R0 = (real_t)pow((double)((nt - 1.0)/(nt + 1.0)), 2.0);
        real_t R = R0 + (1.0 - R0) * (real_t)pow(1.0 - c, 5.0);
        Ray refract_r = Ray(inter.world_intersection, normalize(t));
         
        if (entering) {
            // Russian Roulette sampling
            if (random() < R) {
                return trace_photon(scene, reflect_r, ph, depth + 1, SPECULAR);
            }
            else {
                return trace_photon(scene, refract_r, ph, depth + 1, SPECULAR);
            }
        }
        else {
            // Russian Roulette sampling
            if (random() < R) {
                return trace_photon(scene, reflect_r, ph, depth + 1, SPECULAR);
            }
            else {
                return trace_photon(scene, refract_r, ph, depth + 1, SPECULAR);
            }
        }
    }
}

/**
 * Emits photons to the scene and collect them in photon maps.
 * @param scene The scene to trace.
 * @return void.
 */
void emit_photons(const Scene* scene) {
    int j;
    const SphereLight light = scene->get_lights()[0];
    for (j = 0; j < NUM_PHOTON; j++) {
        Ray r = photon_ray(light);
	// scale the intensity of the photon ray
        Photon ph = Photon(light.color * FACTOR);
	// do photon tracing
        trace_photon(scene, r, ph, 0, 0);
    }   
    
    // balance the maps with kd tree structure
    construct_kdtree(global_illum_map, 0, global_illum_map.size());
    construct_kdtree(caustic_map, 0, caustic_map.size());
}

/**
 * Performs a raytrace on the given pixel on the current scene.
 * The pixel is relative to the bottom-left corner of the image.
 * @param scene The scene to trace.
 * @param x The x-coordinate of the pixel to trace.
 * @param y The y-coordinate of the pixel to trace.
 * @param width The width of the screen in pixels.
 * @param height The height of the screen in pixels.
 * @return The color of that pixel in the final image.
 */
Color3 Raytracer::trace_pixel(const Scene* scene,
			      size_t x,
			      size_t y,
			      size_t width,
			      size_t height,
                              Neighbor *neighbors)
{
    assert(x < width);
    assert(y < height);
    
    real_t dx = real_t(1)/width;
    real_t dy = real_t(1)/height;   
    Vector3 e = scene->camera.get_position(); 
    Color3 res = Color3::Black();
    unsigned int iter;
    for (iter = 0; iter < num_samples; iter++)
    {
        // compute the ray shooting our from the pixel
        real_t i = real_t(2)*(real_t(x) + random())*dx - real_t(1);
        real_t j = real_t(2)*(real_t(y) + random())*dy - real_t(1);
        Vector3 dir = Ray::get_pixel_dir(i, j);
        Ray r = Ray(e, dir);
        
        res += clamp(get_pixel_color(scene, r, 0, neighbors), 0.0, 1.0);

    }
    return res*(real_t(1)/num_samples);
}

/**
 * Raytraces some portion of the scene. Should raytrace for about
 * max_time duration and then return, even if the raytrace is not copmlete.
 * The results should be placed in the given buffer.
 * @param buffer The buffer into which to place the color data. It is
 *  32-bit RGBA (4 bytes per pixel), in row-major order.
 * @param max_time, If non-null, the maximum suggested time this
 *  function raytrace before returning, in seconds. If null, the raytrace
 *  should run to completion.
 * @return true if the raytrace is complete, false if there is more
 *  work to be done.
 */
bool Raytracer::raytrace(unsigned char* buffer, real_t* max_time)
{
    // TODO Add any modifications to this algorithm, if needed.
    
    static const size_t PRINT_INTERVAL = 64;

    // the time in milliseconds that we should stop
    unsigned int end_time = 0;
    bool is_done;

    if (max_time)
    {
        // convert duration to milliseconds
        unsigned int duration = (unsigned int) (*max_time * 1000);
        end_time = SDL_GetTicks() + duration;
    }

    // until time is up, run the raytrace. we render an entire group of
    // rows at once for simplicity and efficiency.
    for (; !max_time || end_time > SDL_GetTicks(); current_row += STEP_SIZE)
    {
        // we're done if we finish the last row
        is_done = current_row >= height;
        // break if we finish
        if (is_done) break;

        int loop_upper = std::min(current_row + STEP_SIZE, height);

        // This tells OpenMP that this loop can be parallelized.
#pragma omp parallel for
        for (int c_row = current_row; c_row < loop_upper; c_row++)
        {
            /*
             * This defines a critical region of code that should be
             * executed sequentially.
             */
#pragma omp critical
            {
                if (c_row % PRINT_INTERVAL == 0)
                    printf("Raytracing (Row %d)\n", c_row);
            }
            // allocate a reusable neighbor list for max heap, instead of
            // allocating new neighbor list for every neighbors lookup          
            Neighbor neighbors[NUM_PHOTON_RADIANCE];
            for (size_t x = 0; x < width; x++)
            {
                // trace a pixel
                Color3 color = trace_pixel(scene, x, c_row, width, height, neighbors);
                // write the result to the buffer, always use 1.0 as the alpha
                color.to_array(&buffer[4 * (c_row * width + x)]);
            }
        }
    }

    if (is_done) printf("Done raytracing!\n");

    return is_done;
}

} /* _462 */
