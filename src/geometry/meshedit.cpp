#include "halfedge.h"

#include <set>
#include <map>
#include <vector>
#include <string>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <spdlog/spdlog.h>

using Eigen::Matrix3f;
using Eigen::Matrix4f;
using Eigen::Vector3f;
using Eigen::Vector4f;
using std::optional;
using std::set;
using std::size_t;
using std::string;
using std::unordered_map;
using std::vector;

HalfedgeMesh::EdgeRecord::EdgeRecord(unordered_map<Vertex*, Matrix4f>& vertex_quadrics, Edge* e)
    : edge(e)
{
    (void)vertex_quadrics;
    optimal_pos = Vector3f(0.0f, 0.0f, 0.0f);
    cost        = 0.0f;
}

bool operator<(const HalfedgeMesh::EdgeRecord& a, const HalfedgeMesh::EdgeRecord& b)
{
    return a.cost < b.cost;
}

optional<Edge*> HalfedgeMesh::flip_edge(Edge* e)
{
    if(e->on_boundary()){
        return std::nullopt;
    }
    // 要用到的半边
    Halfedge* h = e->halfedge;
    Halfedge* h_inv = h->inv;
    Halfedge* h_2_3 = h->next;
    Halfedge* h_3_1 = h_2_3->next;
    Halfedge* h_1_4 = h_inv->next;
    Halfedge* h_4_2 = h_1_4->next;
    // 要用到的顶点
    // v1 and v2 are vertices along the edge
    Vertex* v1 = h->from;
    Vertex* v2 = h_inv->from;
    // v3 and v4 are vertices opposite the edge
    Vertex* v3 = h_3_1->from;
    Vertex* v4 = h_4_2->from;
    // 要用到的面片
    Face* f1 = h->face;
    Face* f2 = h_inv->face;
    // // 重新连接各基本元素
    // h->next = h_3_1;
    // h->prev = h_1_4;
    // h->from = v4;

    // // 其余部分请自己完成
    // h_inv->next = h_4_2;
    // h_inv->prev = h_2_3;
    // h_inv->from = v3;

    // h_3_1->next = h_1_4;
    // h_3_1->prev = h;

    // h_1_4->next = h;
    // h_1_4->prev = h_3_1;
    // h_1_4->face = f1;

    // h_2_3->next = h_inv;
    // h_2_3->prev = h_4_2;
    // h_2_3->face = f2;

    // h_4_2->next = h_2_3;
    // h_4_2->prev = h_inv;

    h->set_neighbors(h_3_1, h_1_4, h->inv, v4, e, f1);
    h_inv->set_neighbors(h_4_2, h_2_3, h_inv->inv, v3, e, f2);
    h_3_1->set_neighbors(h_1_4, h, h_3_1->inv, v3, h_3_1->edge, f1);
    h_1_4->set_neighbors(h, h_3_1, h_1_4->inv, v1, h_1_4->edge, f1);
    h_2_3->set_neighbors(h_inv, h_4_2, h_2_3->inv, v2, h_2_3->edge, f2);
    h_4_2->set_neighbors(h_2_3, h_inv, h_4_2->inv, v4, h_4_2->edge, f2);

    f1->halfedge = h;
    f2->halfedge = h_inv;

    return e;
}

optional<Vertex*> HalfedgeMesh::split_edge(Edge* e)
{
    // 要用到的半边
    Halfedge* h = e->halfedge;
    Halfedge* h_inv = h->inv;
    Halfedge* h_2_3 = h->next;
    Halfedge* h_3_1 = h_2_3->next;
    Halfedge* h_1_4 = h_inv->next;
    Halfedge* h_4_2 = h_1_4->next;
    // 要用到的顶点
    // v1 and v2 are vertices along the edge
    Vertex* v1 = h->from;
    Vertex* v2 = h_inv->from;
    // v3 and v4 are vertices opposite the edge
    Vertex* v3 = h_3_1->from;
    Vertex* v4 = h_4_2->from;
    // 要用到的面片
    Face* f123 = h->face;
    Face* f214 = h_inv->face;
    //要创建的元素
    Vertex* v_new = new_vertex();
    Edge* e2      = new_edge();
    Edge* e3      = new_edge();
    Edge* e4      = new_edge();
    Face* fn14 = new_face();
    Face* fn23 = new_face();
    Halfedge* h_n_3 = new_halfedge();
    Halfedge* h_3_n = new_halfedge();
    Halfedge* h_n_2 = new_halfedge();
    Halfedge* h_n_4 = new_halfedge();
    Halfedge* h_4_n = new_halfedge();
    Halfedge* h_n_1 = new_halfedge();
    // 重新连接各基本元素
    v_new->halfedge   = h_n_1;
    v_new->pos        = e->center();

    h_n_3->set_neighbors(h_3_1, h, h_3_n, v_new, e3, f123);
    h_3_n->set_neighbors(h_n_2, h_2_3, h_n_3, v3, e3, fn23);
    h_n_2->set_neighbors(h_2_3, h_3_n, h_inv, v_new, e2, fn23);
    h_inv->set_neighbors(h_n_4, h_4_2, h_n_2, v2, e2, f214);
    h_n_4->set_neighbors(h_4_2, h_inv, h_4_n, v_new, e4, f214);
    h_4_n->set_neighbors(h_n_1, h_1_4, h_n_4, v4, e4, fn14);
    h_n_1->set_neighbors(h_1_4, h_4_n, h, v_new, e, fn14);
    h->set_neighbors(h_n_3, h_3_1, h_n_1, v1, e, f123);

    e->halfedge = h;
    e2->halfedge = h_inv;
    e3->halfedge = h_3_n;
    e4->halfedge = h_4_n;
    e3->is_new     = true;
    e4->is_new     = true;
    f123->halfedge = h_3_1;
    f214->halfedge = h_4_2;
    fn23->halfedge = h_2_3;
    fn14->halfedge = h_1_4;

    h_3_1->set_neighbors(h, h_n_3, h_3_1->inv, v3, h_3_1->edge, f123);
    h_1_4->set_neighbors(h_4_n, h_n_1, h_1_4->inv, v1, h_1_4->edge, fn14);
    h_4_2->set_neighbors(h_inv, h_n_4, h_4_2->inv, v4, h_4_2->edge, f214);
    h_2_3->set_neighbors(h_3_n, h_n_2, h_2_3->inv, v2, h_2_3->edge, fn23);

    return v_new;
}

optional<Vertex*> HalfedgeMesh::collapse_edge(Edge* e)
{
    if(!e->on_boundary()){
        // 处理非边界情况
        // 要用到的半边
        Halfedge* h_1_2 = e->halfedge;
        Halfedge* h_2_1 = h_1_2->inv;
        Halfedge* h_2_3 = h_1_2->next;
        Halfedge* h_3_1 = h_2_3->next;
        Halfedge* h_1_4 = h_2_1->next;
        Halfedge* h_4_2 = h_1_4->next;
        
        Halfedge* h1    = h_3_1->inv;
        Halfedge* h_inv_1 = h_2_3->inv;
        Halfedge* h2    = h_4_2->inv;
        Halfedge* h_inv_2 = h_1_4->inv;
        // 要用到的顶点
        // v1 and v2 are vertices along the edge
        Vertex* v1 = h_1_2->from;
        Vertex* v2 = h_2_1->from;
        // v3 and v4 are vertices opposite the edge
        Vertex* v3 = h_3_1->from;
        Vertex* v4 = h_4_2->from;
        // 要用到的边
        Edge* e1 = h1->edge;
        Edge* e_inv_1 = h_inv_1->edge;
        Edge* e2 = h2->edge;
        Edge* e_inv_2 = h_inv_2->edge;
        // 要用到的面片
        Face* f1 = h_1_2->face;
        Face* f2 = h_2_1->face;
        // 要创建的元素
        Vertex* v_new = new_vertex();
        v_new->pos    = e->center();

        // 将原本连接到 v1 和 v2 的半边连接到v_new
        Halfedge* h = v1->halfedge;
        do {
            h->from = v_new;
            h = h->inv->next;
        } while (h != v1->halfedge);

        h = v2->halfedge;
        do {
            h->from = v_new;
            h = h->inv->next;
        } while (h != v2->halfedge);

        h1->set_neighbors(h1->next, h1->prev, h_inv_1, v_new, e1, h1->face);
        h_inv_1->set_neighbors(h_inv_1->next, h_inv_1->prev, h1, v3, e1, h_inv_1->face);
        h2->set_neighbors(h2->next, h2->prev, h_inv_2, v_new, e2, h2->face);
        h_inv_2->set_neighbors(h_inv_2->next, h_inv_2->prev, h2, v4, e2, h_inv_2->face);
        e1->halfedge = h1;
        e2->halfedge = h2;
        v_new->halfedge = h1;
        v3->halfedge    = h_inv_1;
        v4->halfedge    = h_inv_2;
        
        // 删除元素
        erase(f1);
        erase(f2);
        erase(v1);
        erase(v2);
        erase(e);
        erase(e_inv_1);
        erase(e_inv_2);
        erase(h_1_2);
        erase(h_2_1);
        erase(h_3_1);
        erase(h_4_2);
        erase(h_1_4);
        erase(h_2_3);
        return v_new;
    }
    else{
        // 处理边界情况
        // 要用到的半边
        Halfedge* h_1_2 = e->halfedge;
        Halfedge* h_2_1 = h_1_2->inv;

        // 确定哪个半边是边界半边
        Halfedge* boundary_h = h_1_2->is_boundary() ? h_1_2 : h_2_1;
        Halfedge* inner_h = boundary_h->inv;

        // 邻接的半边
        Halfedge* h_next = inner_h->next;
        Halfedge* h_prev = inner_h->prev;
        Halfedge* h1 = h_prev->inv;
        Halfedge* h_inv_1 = h_next->inv;

        // 要用到的顶点
        Vertex* v1 = inner_h->from;
        Vertex* v2 = boundary_h->from;
        Vertex* v3 = h_prev->from;

        // 要用到的边
        Edge* e1 = h1->edge;
        Edge* e_inv_1 = h_inv_1->edge;

        // 要用到的面片
        Face* f1 = inner_h->face;

        // 要创建的新顶点
        Vertex* v_new = new_vertex();
        v_new->pos = e->center();
        v_new->halfedge = boundary_h;

        // 更新边界半边的起点为 v_new
        boundary_h->from = v_new;

        // 更新连接到 v1 和 v2 的半边，使其指向新顶点 v_new
        Halfedge* h = v1->halfedge;
        do {
            h->from = v_new;
            h = h->inv->next;
        } while (h != v1->halfedge);

        h = v2->halfedge;
        do {
            h->from = v_new;
            h = h->inv->next;
        } while (h != v2->halfedge);

        // 更新内部半边的关系
        h1->set_neighbors(h1->next, h1->prev, h_inv_1, v_new, e1, h1->face);
        h_inv_1->set_neighbors(h_inv_1->next, h_inv_1->prev, h1, v3, e1, h_inv_1->face);
        e1->halfedge = h1;
        v3->halfedge = h_inv_1;

        // 删除不需要的元素
        erase(v1);
        erase(v2);
        erase(e);
        erase(e_inv_1);
        erase(f1);
        erase(inner_h);
        erase(h_next);
        erase(h_prev);

        return v_new;
    }
}

void HalfedgeMesh::loop_subdivide()
{
    optional<HalfedgeMeshFailure> check_result = validate();
    if (check_result.has_value()) {
        return;
    }
    logger->info("subdivide object {} (ID: {}) with Loop Subdivision strategy", object.name,
                 object.id);
    logger->info("original mesh: {} vertices, {} faces in total", vertices.size, faces.size);
    // Each vertex and edge of the original mesh can be associated with a vertex
    // in the new (subdivided) mesh.
    // Therefore, our strategy for computing the subdivided vertex locations is to
    // *first* compute the new positions using the connectivity of the original
    // (coarse) mesh. Navigating this mesh will be much easier than navigating
    // the new subdivided (fine) mesh, which has more elements to traverse.
    // We will then assign vertex positions in the new mesh based on the values
    // we computed for the original mesh.
    
    // Compute new positions for all the vertices in the input mesh using
    // the Loop subdivision rule and store them in Vertex::new_pos.
    //    At this point, we also want to mark each vertex as being a vertex of the
    //    original mesh. Use Vertex::is_new for this.
    for (Vertex* v = vertices.head; v != nullptr; v = v->next_node) {
        {
            int n      = int(v->degree());
            float u    = (n == 3) ? (3.0f / 16.0f) : (3.0f / (8.0f * n));
            v->new_pos = (1.0f - n * u) * v->pos + u * n * v->neighborhood_center();
        }
    }
    // Next, compute the subdivided vertex positions associated with edges, and
    // store them in Edge::new_pos.
    for (Edge* e = edges.head; e != nullptr; e = e->next_node) {
        {
            Vertex* v1 = e->halfedge->from;
            Vertex* v2 = e->halfedge->inv->from;
            
            Vertex* v3 = e->halfedge->next->inv->from;
            Vertex* v4 = e->halfedge->inv->next->inv->from;

            e->new_pos = (3.0f / 8.0f) * (v1->pos + v2->pos) + (1.0f / 8.0f) * (v3->pos + v4->pos);
        }
    }
    // Next, we're going to split every edge in the mesh, in any order.
    // We're also going to distinguish subdivided edges that came from splitting
    // an edge in the original mesh from new edges by setting the boolean Edge::is_new.
    // Note that in this loop, we only want to iterate over edges of the original mesh.
    // Otherwise, we'll end up splitting edges that we just split (and the
    // loop will never end!)
    // I use a vector to store iterators of original because there are three kinds of
    // edges: original edges, edges split from original edges and newly added edges.
    // The newly added edges are marked with Edge::is_new property, so there is not
    // any other property to mark the edges I just split.
    std::vector<Edge*> originalEdges;
    for (Edge* e = edges.head; e != nullptr; e = e->next_node) {
        originalEdges.push_back(e);
    }
    for (Edge* e : originalEdges) {
        auto optional_vertex = split_edge(e);
        if (optional_vertex) {
            Vertex* v = *optional_vertex;
            v->pos = e->new_pos;
            v->is_new = true;
        }
    }
    // Now flip any new edge that connects an old and new vertex.
    for (Edge* e = edges.head; e != nullptr; e = e->next_node) {
        Halfedge* h = e->halfedge;
        if (e->is_new && (h->from->is_new ^ h->inv->from->is_new)){
            flip_edge(e);
        }
    }
    // Finally, copy new vertex positions into the Vertex::pos.
    for (Vertex* v = vertices.head; v != nullptr; v = v->next_node) {
        if (!v->is_new){
            v->pos = v->new_pos;
        }
        v->is_new = false;
    }
    for (Edge* e = edges.head; e != nullptr; e = e->next_node){
        e->is_new = false;
    }
    // Once we have successfully subdivided the mesh, set global_inconsistent
    // to true to trigger synchronization with GL::Mesh.
    global_inconsistent = true;
    logger->info("subdivided mesh: {} vertices, {} faces in total", vertices.size, faces.size);
    logger->info("Loop Subdivision done");
    logger->info("");
    validate();
}

void HalfedgeMesh::simplify()
{
    optional<HalfedgeMeshFailure> check_result = validate();
    if (check_result.has_value()) {
        return;
    }
    logger->info("simplify object {} (ID: {})", object.name, object.id);
    logger->info("original mesh: {} vertices, {} faces", vertices.size, faces.size);
    unordered_map<Vertex*, Matrix4f> vertex_quadrics;
    unordered_map<Face*, Matrix4f> face_quadrics;
    unordered_map<Edge*, EdgeRecord> edge_records;
    set<EdgeRecord> edge_queue;

    // Compute initial quadrics for each face by simply writing the plane equation
    // for the face in homogeneous coordinates. These quadrics should be stored
    // in face_quadrics

    // -> Compute an initial quadric for each vertex as the sum of the quadrics
    //    associated with the incident faces, storing it in vertex_quadrics

    // -> Build a priority queue of edges according to their quadric error cost,
    //    i.e., by building an Edge_Record for each edge and sticking it in the
    //    queue. You may want to use the above PQueue<Edge_Record> for this.

    // -> Until we reach the target edge budget, collapse the best edge. Remember
    //    to remove from the queue any edge that touches the collapsing edge
    //    BEFORE it gets collapsed, and add back into the queue any edge touching
    //    the collapsed vertex AFTER it's been collapsed. Also remember to assign
    //    a quadric to the collapsed vertex, and to pop the collapsed edge off the
    //    top of the queue.

    logger->info("simplified mesh: {} vertices, {} faces", vertices.size, faces.size);
    logger->info("simplification done\n");
    global_inconsistent = true;
    validate();
}

void HalfedgeMesh::isotropic_remesh()
{
    optional<HalfedgeMeshFailure> check_result = validate();
    if (check_result.has_value()) {
        return;
    }
    logger->info("remesh the object {} (ID: {}) with strategy Isotropic Remeshing", object.name,
                 object.id);
    logger->info("original mesh: {} vertices, {} faces", vertices.size, faces.size);
    // Compute the mean edge length.

    // Repeat the four main steps for 5 or 6 iterations
    // -> Split edges much longer than the target length (being careful about
    //    how the loop is written!)
    // -> Collapse edges much shorter than the target length.  Here we need to
    //    be EXTRA careful about advancing the loop, because many edges may have
    //    been destroyed by a collapse (which ones?)
    // -> Now flip each edge if it improves vertex degree
    // -> Finally, apply some tangential smoothing to the vertex positions
    static const size_t iteration_limit = 5;
    set<Edge*> selected_edges;
    for (size_t i = 0; i != iteration_limit; ++i) {
        // Split long edges.

        // Collapse short edges.

        // Flip edges.

        // Vertex averaging.
    }
    logger->info("remeshed mesh: {} vertices, {} faces\n", vertices.size, faces.size);
    global_inconsistent = true;
    validate();
}
