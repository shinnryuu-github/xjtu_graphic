#include "ray.h"

#include <cmath>
#include <array>

#include <Eigen/Dense>
#include <spdlog/spdlog.h>

#include "../utils/math.hpp"

using Eigen::Matrix3f;
using Eigen::Matrix4f;
using Eigen::Vector2f;
using Eigen::Vector3f;
using Eigen::Vector4f;
using std::numeric_limits;
using std::optional;
using std::size_t;

constexpr float infinity = 1e5f;
constexpr float eps      = 1e-5f;

Intersection::Intersection() : t(numeric_limits<float>::infinity()), face_index(0)
{
}

Ray generate_ray(int width, int height, int x, int y, Camera& camera, float depth)
{
    // these lines below are just for compiling and can be deleted
    (void)width;
    (void)height;
    (void)x;
    (void)y;
    (void)depth;
    // these lines above are just for compiling and can be deleted


    // The ratio between the specified plane (width x height)'s depth and the image plane's depth.
    
    // Transfer the view-space position to world space.
    Vector3f world_pos;
    return {camera.position, (world_pos - camera.position).normalized()};
}

optional<Intersection> ray_triangle_intersect(const Ray& ray, const GL::Mesh& mesh, size_t index)
{
    // these lines below are just for compiling and can be deleted
    (void)ray;
    (void)mesh;
    (void)index;
    // these lines above are just for compiling and can be deleted
    Intersection result;
    
    if (result.t - infinity < -eps) {
        return result;
    } else {
        return std::nullopt;
    }
}

optional<Intersection> naive_intersect(const Ray& ray, const GL::Mesh& mesh, const Matrix4f model)
{
    Intersection result;
    result.t = infinity;
    for (size_t i = 0; i < mesh.faces.count(); ++i) {
        // Vertex a, b and c are assumed to be in counterclockwise order.
        // Construct matrix A = [d, a - b, a - c] and solve Ax = (a - origin)
        // Matrix A is not invertible, indicating the ray is parallel with the triangle.
        // Test if alpha, beta and gamma are all between 0 and 1.

        // 获取三角形的顶点
        const auto& face = mesh.face(i);
        Eigen::Vector4f a = model * Eigen::Vector4f(mesh.vertex(face[0])(0), mesh.vertex(face[0])(1), mesh.vertex(face[0])(2), 1.0f);
        Eigen::Vector4f b = model * Eigen::Vector4f(mesh.vertex(face[1])(0), mesh.vertex(face[1])(1), mesh.vertex(face[1])(2), 1.0f);
        Eigen::Vector4f c = model * Eigen::Vector4f(mesh.vertex(face[2])(0), mesh.vertex(face[2])(1), mesh.vertex(face[2])(2), 1.0f);
        // 提取 3D 坐标
        Eigen::Vector3f a3 = a.head<3>();
        Eigen::Vector3f b3 = b.head<3>();
        Eigen::Vector3f c3 = c.head<3>();

        // 构造向量
        Eigen::Vector3f ab     = b3 - a3;
        Eigen::Vector3f ac     = c3 - a3;
        Eigen::Vector3f d      = ray.direction;
        Eigen::Vector3f origin = ray.origin;

        // 计算行列式
        Eigen::Matrix3f A;
        A.col(0) = -d;
        A.col(1) = ab;
        A.col(2) = ac;

        // 计算行列式
        float det = A.determinant();
        if (std::fabs(det) < eps) {
            // 行列式接近于零，表示光线与三角形平行
            continue;
        }

        // 解线性方程组
        Eigen::Vector3f AO = origin - a3;
        Eigen::Vector3f x = A.inverse() * AO;

        float t = x[0];
        float alpha = x[1];
        float beta = x[2];
        float gamma = 1.0f - alpha - beta;

        if (t < eps) {
            continue;
        }

        // 检查 alpha, beta, gamma 是否都在 [0, 1] 之间
        if (alpha >= -eps && alpha <= 1 + eps &&
            beta >= -eps && beta <= 1 + eps &&
            gamma >= -eps && gamma <= 1 + eps) {
            if (t < result.t) {
                result.t = t;
                result.face_index = i;
                result.barycentric_coord = Eigen::Vector3f(alpha, beta, gamma);
                result.normal = ab.cross(ac).normalized();
            }
        }
    }
    // Ensure result.t is strictly less than the constant `infinity`.
    if (result.t - infinity < -eps) {
        return result;
    }
    return std::nullopt;
}
