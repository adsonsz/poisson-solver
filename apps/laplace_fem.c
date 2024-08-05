#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

// Bijection between N and N^2.
struct Pair {
    int a;
    int b;
};

struct Vector {
    double x;
    double y;
};

// Very inneficient Map data structure.
// Perhaps I shall change later on.
struct Map {
    int* key;
    int* value;

    int size;
    int capacity;
};

// Very inneficient set data structure.
// Perhaps I shall change later on.
struct Set {
    int* data;
    int size;
    int capacity;
};

// Triangular mesh.
// Each vertex is of the form (x, y), 2D mesh.
// Each face is a triangle, an index of three vertices.
struct Mesh {
    // Vertex information of the mesh.
    double* vx;
    double* vy;
    int amount_vertices;
    int capacity_vertices;

    // Face information of the mesh. Mesh is triangular.
    // Eg.: f1 (v1, v2, v3)
    int* af;
    int* bf;
    int* cf;
    int amount_faces;
    int capacity_faces;

    // FV data structure.
    // Eg.: v1 (f1, f4, f8, ..)
    int** fv;
    int* amount_fv;
    int* capacity_fv;
};

// Finite element method struct.
struct FEM {
    struct Mesh* mesh;
    double* vertex_data;
    bool* mark_boundary;
};

// Determinant function.
double det(double x1, double y1, double x2, double y2);

// Copy an array into another.
void copy_array(int size, double* from, double* to);
void zero_array(int size, double* arr);

// Pair-number bijection.
int integer_from_pair(struct Pair pair);
int integer_from_pair_direct(int a, int b);
struct Pair pair_from_integer(int num);

// Set struct.
struct Set* set_create(int sz);
void set_free(struct Set* map);
void set_add(struct Set* set, int element);
bool set_contains(struct Set* set, int element);

// Mapping struct.
struct Map* map_create(int sz);
void map_free(struct Map* map);
void map_add_pair(struct Map* map, int key, int value);
bool map_contains(struct Map* map, int key);
int map_get_value(struct Map* map, int key);

// Mesh construction, allocation, reallocation, and destruction.
struct Mesh* mesh_create(int size);
void mesh_free(struct Mesh* mesh);
void mesh_realloc_vertices(struct Mesh* mesh);
void mesh_realloc_faces(struct Mesh* mesh);
void mesh_realloc_fv(struct Mesh* mesh, int index);

// Mesh assembly.
int mesh_add_vertex(struct Mesh* mesh, double x, double y);
int mesh_add_face(struct Mesh* mesh, int a, int b, int c);

// Get mesh information.
double mesh_get_face_area(struct Mesh* mesh, int face);

// Mesh traversal
struct Pair mesh_get_adjacent_faces_at_edge(struct Mesh* mesh, int v1, int v2);
struct Set* mesh_get_connected_vertices_at_vertex(struct Mesh* mesh, int vertex);

// Subdivision algorithm.
struct Mesh* mesh_simple_subdivide(struct Mesh* mesh);

// Save mesh.
void mesh_save(struct Mesh* mesh, char* filename);

// Finite element methods.
struct FEM* fem_create(int iter);
void fem_free(struct FEM* fem);
struct Vector fem_gradient_component(struct FEM* fem, int vertex, int face);
double fem_diagonal_element(struct FEM* fem, int vertex);
double fem_non_diagonal_adjacent(struct FEM* fem, int vertex, int adj);
void fem_jacobi_iteration(struct FEM* fem, double* aux);



// ----------------------------------------------------------------------------------------

double det(double x1, double y1, double x2, double y2) {
    return x1 * y2 - x2 * y1;
}

void copy_array(int size, double* from, double* to) {
    for (int i = 0; i < size; ++i) to[i] = from[i];
}

void zero_array(int size, double* arr) {
    for (int i = 0; i < size; ++i) arr[i] = 0.0;
}


struct Set* set_create(int sz) {
    struct Set* set = malloc(sizeof(struct Set));
    set->data = malloc(sz * sizeof(int));
    set->capacity = sz;
    set->size = 0;

    return set;
}

void set_free(struct Set* set) {
    free(set->data);
    free(set);
}

void set_add(struct Set* set, int element) {
    if (set_contains(set, element) == true) return;

    if (set->size >= set->capacity) {
        set->capacity *= 2;
        set->data = realloc(set->data, set->capacity * sizeof(int));
    }

    set->data[set->size++] = element;
}

bool set_contains(struct Set* set, int element) {
    for (int i = 0; i < set->size; ++i) {
        if (set->data[i] == element) return true;
    }

    return false;
}

struct Map* map_create(int sz) {
    struct Map* map = malloc(sizeof(struct Map));
    map->key = malloc(sz * sizeof(int));
    map->value = malloc(sz * sizeof(int));
    map->capacity = sz;
    map->size = 0;

    return map;
}

void map_free(struct Map* map) {
    free(map->key);
    free(map->value);
    free(map);
}

void map_add_pair(struct Map* map, int key, int value) {
    if (map_contains(map, key) == false) {
        if (map->size >= map->capacity) {
            map->capacity *= 2;
            map->key = realloc(map->key, map->capacity * sizeof(int));
            map->value = realloc(map->value, map->capacity * sizeof(int));
        }

        map->key[map->size] = key;
        map->value[map->size] = value;
        map->size += 1;
    }
}

bool map_contains(struct Map* map, int key) {
    for (int i = 0; i < map->size; ++i) {
        if (map->key[i] == key) return true;
    }

    return false;
}

int map_get_value(struct Map* map, int key) {
    for (int i = 0; i < map->size; ++i) {
        if (map->key[i] == key) return map->value[i];
    }

    return -1;
}

int integer_from_pair(struct Pair pair) {
    return integer_from_pair_direct(pair.a, pair.b);
}

int integer_from_pair_direct(int a, int b) {
    int diag = a + b;
    return diag * (diag + 1) / 2 + b;
}

struct Pair pair_from_integer(int num) {
    float fnum = (float)num;
    float ca = floor((sqrt(8*fnum + 1) - 1) / 2);
    int w = (int)ca;
    int t = (w*w + w) / 2;
    int y = num - t;
    int x = w - y;

    struct Pair pair;
    pair.a = x;
    pair.b = y;

    return pair;
}

struct Mesh* mesh_create(int size) {
    // Allocate mesh.
    struct Mesh* mesh = malloc(sizeof(struct Mesh));

    // Allocate vertices.
    mesh->amount_vertices = 0;
    mesh->capacity_vertices = size;
    mesh->vx = malloc(size * sizeof(double));
    mesh->vy = malloc(size * sizeof(double));

    // Allocate faces.
    mesh->amount_faces = 0;
    mesh->capacity_faces = size;
    mesh->af = malloc(size * sizeof(int));
    mesh->bf = malloc(size * sizeof(int));
    mesh->cf = malloc(size * sizeof(int));

    // Allocate FV data structure.
    mesh->capacity_fv = malloc(size * sizeof(int));
    mesh->amount_fv= malloc(size * sizeof(int));
    mesh->fv = malloc(size * sizeof(int*));
    for (int i = 0; i < mesh->capacity_vertices; ++i) {
        mesh->amount_fv[i] = 0;
        mesh->capacity_fv[i] = 1;
        mesh->fv[i] = malloc(mesh->capacity_fv[i] * sizeof(int));
    }

    return mesh;
}

void mesh_realloc_vertices(struct Mesh* mesh) {
    int prev = mesh->capacity_vertices;
    mesh->capacity_vertices *= 2;
    mesh->vx = realloc(mesh->vx, mesh->capacity_vertices * sizeof(double));
    mesh->vy = realloc(mesh->vy, mesh->capacity_vertices * sizeof(double));

    mesh->fv = realloc(mesh->fv, mesh->capacity_vertices * sizeof(int*));
    mesh->amount_fv = realloc(mesh->amount_fv, mesh->capacity_vertices * sizeof(int));
    mesh->capacity_fv = realloc(mesh->capacity_fv, mesh->capacity_vertices * sizeof(int));

    for (int i = prev; i < mesh->capacity_vertices; ++i) {
        mesh->amount_fv[i] = 0;
        mesh->capacity_fv[i] = 1;
        mesh->fv[i] = malloc(mesh->capacity_fv[i] * sizeof(int));
    }
}

void mesh_realloc_faces(struct Mesh* mesh) {
    mesh->capacity_faces *= 2;
    mesh->af = realloc(mesh->af, mesh->capacity_faces * sizeof(int));
    mesh->bf = realloc(mesh->bf, mesh->capacity_faces * sizeof(int));
    mesh->cf = realloc(mesh->cf, mesh->capacity_faces * sizeof(int));
}

void mesh_realloc_fv(struct Mesh* mesh, int index) {
    mesh->capacity_fv[index] *= 2;
    mesh->fv[index] = realloc(mesh->fv[index], mesh->capacity_fv[index] * sizeof(int));
}

void mesh_free(struct Mesh* mesh) {
    free(mesh->vx);
    free(mesh->vy);
    free(mesh->af);
    free(mesh->bf);
    free(mesh->cf);

    free(mesh->capacity_fv);
    free(mesh->amount_fv);
    for (int i = 0; i < mesh->capacity_vertices; ++i) {
        free(mesh->fv[i]);
    }

    free(mesh->fv);
    free(mesh);
}

int mesh_add_vertex(struct Mesh* mesh, double x, double y) {
    int sz = mesh->amount_vertices;
    if (sz >=  mesh->capacity_vertices) {
        mesh_realloc_vertices(mesh);
    }

    mesh->vx[sz] = x;
    mesh->vy[sz] = y;
    mesh->amount_vertices += 1;

    return sz;
}

int mesh_add_face(struct Mesh* mesh, int v1, int v2, int v3) {
    int face_index = mesh->amount_faces;
    if (mesh->amount_faces >=  mesh->capacity_faces) {
        mesh_realloc_faces(mesh);
    }

    mesh->af[face_index] = v1;
    mesh->bf[face_index] = v2;
    mesh->cf[face_index] = v3;
    mesh->amount_faces += 1;

    // Register the FV data structure.
    if (mesh->amount_fv[v1] >= mesh->capacity_fv[v1]) mesh_realloc_fv(mesh, v1);
    if (mesh->amount_fv[v2] >= mesh->capacity_fv[v2]) mesh_realloc_fv(mesh, v2);
    if (mesh->amount_fv[v3] >= mesh->capacity_fv[v3]) mesh_realloc_fv(mesh, v3);

    // Recall.: v1 (f1, f4, f8, ..)
    mesh->fv[v1][mesh->amount_fv[v1]++] = face_index;
    mesh->fv[v2][mesh->amount_fv[v2]++] = face_index;
    mesh->fv[v3][mesh->amount_fv[v3]++] = face_index;

    return face_index;
}

double mesh_get_face_area(struct Mesh* mesh, int face) {
    int v1 = mesh->af[face];
    int v2 = mesh->bf[face];
    int v3 = mesh->cf[face];

    double vx1 = mesh->vx[v1];
    double vy1 = mesh->vy[v1];
    double vx2 = mesh->vx[v2];
    double vy2 = mesh->vy[v2];
    double vx3 = mesh->vx[v3];
    double vy3 = mesh->vy[v3];

    double signed_area = 0.5 * (
        det(vx3, vy3, vx1, vy1)
        + det(vx1, vy1, vx2, vy2)
        + det(vx2, vy2, vx3, vy3)
    );

    if (signed_area < 0) return -signed_area;
    else return signed_area;
}

struct Pair mesh_get_adjacent_faces_at_edge(struct Mesh* mesh, int v1, int v2) {
    struct Pair pair;
    pair.a = -1;
    pair.b = -1;

    for (int i = 0; i < mesh->amount_fv[v1]; ++i) {
        if (pair.a != -1  &&  pair.b != -1)  break;
        for (int j = 0; j < mesh->amount_fv[v2]; ++j) {
            if (pair.a != -1  &&  pair.b != -1)  break;

            int face1 = mesh->fv[v1][i];
            int face2 = mesh->fv[v2][j];

            if (face1 == face2) {
                if (pair.a == -1) {
                    pair.a = face1;
                } else {
                    if (pair.a == face1) continue;
                    else pair.b = face1;
                }
            }
        }
    }

    return pair;
}

struct Set* mesh_get_connected_vertices_at_vertex(struct Mesh* mesh, int vertex) {
    struct Set* set = set_create(5);

    for (int i = 0; i < mesh->amount_fv[vertex]; ++i) {
        int face = mesh->fv[vertex][i];
        if (mesh->af[face] == vertex) {
            set_add(set, mesh->bf[face]);
            set_add(set, mesh->cf[face]);
        }

        if (mesh->bf[face] == vertex) {
            set_add(set, mesh->cf[face]);
            set_add(set, mesh->af[face]);
        }

        if (mesh->cf[face] == vertex) {
            set_add(set, mesh->af[face]);
            set_add(set, mesh->bf[face]);
        }
    }

    return set;
}

struct Mesh* mesh_simple_subdivide(struct Mesh* mesh) {

    int estimated_amount_vertices = 3 * mesh->amount_vertices + 5;
    struct Mesh* result = mesh_create(estimated_amount_vertices);

    // Copy vertices into new mesh.
    for (int i = 0; i < mesh->amount_vertices; ++i) {
        mesh_add_vertex(result, mesh->vx[i], mesh->vy[i]);
    }

    // Save information. [Edge (v1, v2), midpoint].
    struct Map* edge_midvertex = map_create(3 * mesh->amount_vertices);
    for (int face = 0; face < mesh->amount_faces; ++face) {
        int v1 = mesh->af[face];
        int v2 = mesh->bf[face];
        int v3 = mesh->cf[face];
    
        int edge12 = integer_from_pair_direct(v1, v2);
        int edge21 = integer_from_pair_direct(v2, v1);
        if (map_contains(edge_midvertex, edge21) == false) {
            double midx = 0.5 * (mesh->vx[v1] + mesh->vx[v2]);
            double midy = 0.5 * (mesh->vy[v1] + mesh->vy[v2]);
            int midvertex = mesh_add_vertex(result, midx, midy);
            map_add_pair(edge_midvertex, edge12, midvertex);
            map_add_pair(edge_midvertex, edge21, midvertex);
        }

        int edge23 = integer_from_pair_direct(v2, v3);
        int edge32 = integer_from_pair_direct(v3, v2);
        if (map_contains(edge_midvertex, edge23) == false) {
            double midx = 0.5 * (mesh->vx[v2] + mesh->vx[v3]);
            double midy = 0.5 * (mesh->vy[v2] + mesh->vy[v3]);
            int midvertex = mesh_add_vertex(result, midx, midy);
            map_add_pair(edge_midvertex, edge23, midvertex);
            map_add_pair(edge_midvertex, edge32, midvertex);
        }

        int edge31 = integer_from_pair_direct(v3, v1);
        int edge13 = integer_from_pair_direct(v1, v3);
        if (map_contains(edge_midvertex, edge31) == false) {
            double midx = 0.5 * (mesh->vx[v3] + mesh->vx[v1]);
            double midy = 0.5 * (mesh->vy[v3] + mesh->vy[v1]);
            int midvertex = mesh_add_vertex(result, midx, midy);
            map_add_pair(edge_midvertex, edge31, midvertex);
            map_add_pair(edge_midvertex, edge13, midvertex);
        }
    }

    // Make the subdivision.
    for (int face = 0; face < mesh->amount_faces; ++face) {
        int v1 = mesh->af[face];
        int v2 = mesh->bf[face];
        int v3 = mesh->cf[face];
        
        int m12 = map_get_value(edge_midvertex, integer_from_pair_direct(v1, v2));
        int m23 = map_get_value(edge_midvertex, integer_from_pair_direct(v2, v3));
        int m31 = map_get_value(edge_midvertex, integer_from_pair_direct(v3, v1));


        mesh_add_face(result, m12, v2, m23);
        mesh_add_face(result, m23, v3, m31);
        mesh_add_face(result, m31, v1, m12);
        mesh_add_face(result, m12, m23, m31);
    }


    // Free used memory and return.
    map_free(edge_midvertex);
    return result;
}

void mesh_save(struct Mesh* mesh, char* filename) {
    FILE* fp = fopen(filename, "w");
    if (fp == NULL) exit(1);

    for (int v = 0; v < mesh->amount_vertices; ++v) {
        fprintf(fp, "v %f %f %f\n", mesh->vx[v], mesh->vy[v], 0.0);
    }

    fprintf(fp, "\n");

    for (int f = 0; f < mesh->amount_faces; ++f) {
        fprintf(fp, "f %d %d %d\n", 1+mesh->af[f], 1+mesh->bf[f], 1+mesh->cf[f]);
    }

    fclose(fp);
}


struct FEM* fem_create(int iter) {
    struct Mesh* mesh = mesh_create(5);
    int v1 = mesh_add_vertex(mesh, 0.0, 0.0);
    int v2 = mesh_add_vertex(mesh, 1.0, 0.0);
    int v3 = mesh_add_vertex(mesh, 1.0, 1.0);
    int v4 = mesh_add_vertex(mesh, 0.0, 1.0);
    int vc = mesh_add_vertex(mesh, 0.5, 0.5);


    mesh_add_face(mesh, v1, vc, v2);
    mesh_add_face(mesh, v2, vc, v3);
    mesh_add_face(mesh, v3, vc, v4);
    mesh_add_face(mesh, v4, vc, v1);

    // Subdivide the mesh an 'iter' amount of iteration.
    struct Mesh* better_mesh;
    for (int i = 1; i < iter; ++i) {
        better_mesh = mesh_simple_subdivide(mesh);
        mesh_free(mesh);
        mesh = better_mesh;
    }

    // Set up auxiliary data.
    double* vertex_data = malloc(mesh->amount_vertices * sizeof(double));
    bool* mark_boundary = malloc(mesh->amount_vertices * sizeof(bool));


    // Set up boundary condition.
    for (int i = 0; i < mesh->amount_vertices; ++i) mark_boundary[i] = false;
    for (int i = 0; i < mesh->amount_vertices; ++i) {
        if (mesh->vy[i] < 1e-4) { mark_boundary[i] = true;  vertex_data[i] = 1.0; }
        else if (mesh->vy[i] > 1-1e-4) { mark_boundary[i] = true;  vertex_data[i] = 0.0; }
        else if (mesh->vx[i] < 1e-4) { mark_boundary[i] = true;  vertex_data[i] = 0.0; }
        else if (mesh->vx[i] > 1-1e-4) { mark_boundary[i] = true;  vertex_data[i] = 0.0; }
        else continue;
    }
    
    // Build the FEM structure.
    struct FEM* fem = malloc(sizeof(struct FEM));
    fem->mesh = mesh;
    fem->vertex_data = vertex_data;
    fem->mark_boundary = mark_boundary;

    return fem;
}

void fem_free(struct FEM* fem) {
    mesh_free(fem->mesh);
    free(fem->vertex_data);
    free(fem);
}

struct Vector fem_gradient_component(struct FEM* fem, int vertex, int face) {
    // Get vertices from that face, in the right order.
    int v1, v2, v3;
    if (fem->mesh->af[face] == vertex) {
        v1 = vertex;
        v2 = fem->mesh->bf[face];
        v3 = fem->mesh->cf[face];
    } else if (fem->mesh->bf[face] == vertex) {
        v1 = vertex;
        v2 = fem->mesh->cf[face];
        v3 = fem->mesh->af[face];
    } else {
        v1 = vertex;
        v2 = fem->mesh->af[face];
        v3 = fem->mesh->bf[face];
    }

    // Calculate vertices position.
    double x2 = fem->mesh->vx[v2] - fem->mesh->vx[v1];
    double x3 = fem->mesh->vx[v3] - fem->mesh->vx[v1];
    double y2 = fem->mesh->vy[v2] - fem->mesh->vy[v1];
    double y3 = fem->mesh->vy[v3] - fem->mesh->vy[v1];
    double xdiff = x2 - x3;
    double ydiff = y2 - y3;

    // Calculate gradient components (c1, c2).
    double den = det(x2, y2, x3, y3);
    double c1 = ydiff / den;
    double c2 = -xdiff / den;

    // Return.
    struct Vector gradient;
    gradient.x = c1;
    gradient.y = c2;
    return gradient;
}

double fem_diagonal_element(struct FEM* fem, int vertex) {
    double result = 0.0;
    int* adj_faces = fem->mesh->fv[vertex];
    for (int i = 0; i < fem->mesh->amount_fv[vertex]; ++i) {
        int face = adj_faces[i];
        struct Vector gradient = fem_gradient_component(fem, vertex, face);
        double dot = gradient.x * gradient.x + gradient.y * gradient.y;
        double area = mesh_get_face_area(fem->mesh, face);
        result += dot * area;
    }

    return result;
}

double fem_non_diagonal_adjacent(struct FEM* fem, int vertex, int adj) {
    double result = 0.0;
    struct Pair face = mesh_get_adjacent_faces_at_edge(fem->mesh, vertex, adj);

    if (face.a != -1) {
        struct Vector main_gradient = fem_gradient_component(fem, vertex, face.a);
        struct Vector adjacent_gradient = fem_gradient_component(fem, adj, face.a);
        double dot = main_gradient.x * adjacent_gradient.x + main_gradient.y * adjacent_gradient.y;
        double area = mesh_get_face_area(fem->mesh, face.a);
        result += dot * area;
    }

    if (face.b != -1) {
        struct Vector main_gradient = fem_gradient_component(fem, vertex, face.b);
        struct Vector adjacent_gradient = fem_gradient_component(fem, adj, face.b);
        double dot = main_gradient.x * adjacent_gradient.x + main_gradient.y * adjacent_gradient.y;
        double area = mesh_get_face_area(fem->mesh, face.b);
        result += dot * area;
    }

    return result;
}

double fem_get_matrix_element(struct FEM* fem, int v1, int v2) {
    if (fem->mark_boundary[v1] == false) {
        if (v1 == v2) return fem_diagonal_element(fem, v1);
        else {
            struct Pair pair = mesh_get_adjacent_faces_at_edge(fem->mesh, v1, v2);
            if (pair.a == -1  &&  pair.b == -1) return 0.0;
            else return fem_non_diagonal_adjacent(fem, v1, v2);
        }
    } else {
        if (v1 == v2) return 1.0;
        else return 0.0;
    }
}

void fem_matrix_jacobi_iteration(struct FEM* fem, double* aux) {
    copy_array(fem->mesh->amount_vertices, fem->vertex_data, aux);
    zero_array(fem->mesh->amount_vertices, fem->vertex_data);
    for (int i = 0; i < fem->mesh->amount_vertices; ++i) {
        for (int j = 0; j < fem->mesh->amount_vertices; ++j) {
            if (fem->mark_boundary[i] == false) {
                double diag = fem_diagonal_element(fem, i);
                if (i == j) continue;
                fem->vertex_data[i] -= fem_get_matrix_element(fem, i, j) / diag * aux[j];
            } else {
                fem->vertex_data[i] = aux[i];
            }
        }
    }
}

void fem_jacobi_iteration(struct FEM* fem, double* aux) {
    copy_array(fem->mesh->amount_vertices, fem->vertex_data, aux);
    zero_array(fem->mesh->amount_vertices, fem->vertex_data);
    for (int vertex = 0; vertex < fem->mesh->amount_vertices; ++vertex) {
        if (fem->mark_boundary[vertex] == false) {
            // Compute diagonal term.
            double diag = fem_diagonal_element(fem, vertex);

            // Compute nondiagonal contributions.
            double nondiag = 0.0;
            struct Set* conn = mesh_get_connected_vertices_at_vertex(fem->mesh, vertex);
            for (int i = 0; i < conn->size; ++i) {
                int conn_vertex = conn->data[i];
                if (conn_vertex == vertex) continue;
                nondiag += fem_non_diagonal_adjacent(fem, vertex, conn_vertex) * aux[conn_vertex];
            }

            set_free(conn);
            fem->vertex_data[vertex] -= nondiag / diag;
        } else {
            fem->vertex_data[vertex] = aux[vertex];
        }
    }
}

void fem_save_vertex_data(struct FEM* fem, char* filename) {
    FILE* fp = fopen(filename, "w");

    for (int i = 0; i < fem->mesh->amount_vertices; ++i) {
        fprintf(fp, "%f ", fem->vertex_data[i]);
    }

    fclose(fp);
}

int main() {
    for (int i = 1; i < 7; ++i) {
        struct FEM* fem = fem_create(i);
        int sz = fem->mesh->amount_vertices;
        double* aux = malloc(sz * sizeof(double));
        for (int i = 0; i < sz; ++i) fem_jacobi_iteration(fem, aux);

        printf("[%d] %d vertices\n", i, sz);

        char mesh_filename[128];
        char dat_filename[128];

        sprintf(mesh_filename, "mesh/fem%d.obj", i);
        sprintf(dat_filename, "mesh/fem%d.dat", i);

        mesh_save(fem->mesh, mesh_filename);
        fem_save_vertex_data(fem, dat_filename);

        fem_free(fem);
        free(aux);
    }

    return 0;
}
