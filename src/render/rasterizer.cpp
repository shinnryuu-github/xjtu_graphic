#include <array>
#include <limits>
#include <tuple>
#include <vector>
#include <algorithm>
#include <cmath>
#include <mutex>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <spdlog/spdlog.h>

#include "rasterizer.h"
#include "triangle.h"
#include "../utils/math.hpp"

using Eigen::Matrix4f;
using Eigen::Vector2i;
using Eigen::Vector3f;
using Eigen::Vector4f;
using std::fill;
using std::tuple;

void Rasterizer::worker_thread()
{
    while (true) {
        VertexShaderPayload payload;
        Triangle triangle;
        {
            // printf("vertex_finish = %d\n vertex_shader_output_queue.size = %ld\n",
            // Context::vertex_finish, Context::vertex_shader_output_queue.size());
            if (Context::vertex_finish && Context::vertex_shader_output_queue.empty()) {
                Context::rasterizer_finish = true;
                return;
            }
            if (Context::vertex_shader_output_queue.size() < 3) {
                continue;
            }
            std::unique_lock<std::mutex> lock(Context::vertex_queue_mutex);
            if (Context::vertex_shader_output_queue.size() < 3) {
                continue;
            }
            for (size_t vertex_count = 0; vertex_count < 3; vertex_count++) {
                payload = Context::vertex_shader_output_queue.front();
                Context::vertex_shader_output_queue.pop();
                if (vertex_count == 0) {
                    triangle.world_pos[0]    = payload.world_position;
                    triangle.viewport_pos[0] = payload.viewport_position;
                    triangle.normal[0]       = payload.normal;
                } else if (vertex_count == 1) {
                    triangle.world_pos[1]    = payload.world_position;
                    triangle.viewport_pos[1] = payload.viewport_position;
                    triangle.normal[1]       = payload.normal;
                } else {
                    triangle.world_pos[2]    = payload.world_position;
                    triangle.viewport_pos[2] = payload.viewport_position;
                    triangle.normal[2]       = payload.normal;
                }
            }
        }
        rasterize_triangle(triangle);
    }
}

float sign(Eigen::Vector2f p1, Eigen::Vector2f p2, Eigen::Vector2f p3)
{
    return (p1.x() - p3.x()) * (p2.y() - p3.y()) - (p2.x() - p3.x()) * (p1.y() - p3.y());
}

// 给定坐标(x,y)以及三角形的三个顶点坐标，判断(x,y)是否在三角形的内部
bool Rasterizer::inside_triangle(int x, int y, const Vector4f* vertices)
{
    Eigen::Vector3f v0 = {vertices[0].x(), vertices[0].y(), 1.0f};
    Eigen::Vector3f v1 = {vertices[1].x(), vertices[1].y(), 1.0f};
    Eigen::Vector3f v2 = {vertices[2].x(), vertices[2].y(), 1.0f};

    Eigen::Vector3f p(float(x), float(y), 1.0f);

    Eigen::Vector3f c0 = (v1 - v0).cross(p - v0);
    Eigen::Vector3f c1 = (v2 - v1).cross(p - v1);
    Eigen::Vector3f c2 = (v0 - v2).cross(p - v2);

    return (c0.z() >= 0 && c1.z() >= 0 && c2.z() >= 0) || 
           (c0.z() <= 0 && c1.z() <= 0 && c2.z() <= 0);
}

// 给定坐标(x,y)以及三角形的三个顶点坐标，计算(x,y)对应的重心坐标[alpha, beta, gamma]
tuple<float, float, float> Rasterizer::compute_barycentric_2d(float x, float y, const Vector4f* v)
{
    Eigen::Vector2f v0(v[0].x(), v[0].y());
    Eigen::Vector2f v1(v[1].x(), v[1].y());
    Eigen::Vector2f v2(v[2].x(), v[2].y());

    float area = 0.5f * (-v1.y() * v2.x() + v0.y() * (-v1.x() + v2.x()) + v0.x() * (v1.y() - v2.y()) + v1.x() * v2.y());

    float area0 = 0.5f * (-v1.y() * x + v0.y() * (-v1.x() + x) + v0.x() * (v1.y() - y) + v1.x() * y);
    float area1 = 0.5f * (-y * v2.x() + v1.y() * (-v2.x() + x) + v1.x() * (y - v2.y()) + x * v2.y());
    float area2 = area - area0 - area1;

    float alpha = area0 / area;
    float beta = area1 / area;
    float gamma = area2 / area;

    return {alpha, beta, gamma};
}

// 对顶点的某一属性插值
Vector3f Rasterizer::interpolate(float alpha, float beta, float gamma, const Eigen::Vector3f& vert1,
                                 const Eigen::Vector3f& vert2, const Eigen::Vector3f& vert3,
                                 const Eigen::Vector3f& weight, const float& Z)
{
    Vector3f interpolated_res;
    for (int i = 0; i < 3; i++) {
        interpolated_res[i] = alpha * vert1[i] / weight[0] + beta * vert2[i] / weight[1] +
                              gamma * vert3[i] / weight[2];
    }
    interpolated_res *= Z;
    return interpolated_res;
}

// 对当前三角形进行光栅化
void Rasterizer::rasterize_triangle(Triangle& t)
{
    // if current pixel is in current triange:
    // 1. interpolate depth(use projection correction algorithm)
    // 2. interpolate vertex positon & normal(use function:interpolate())
    // 3. push primitive into fragment queue

    // 计算三角形的边界框
    int min_x = std::min({static_cast<int>(t.viewport_pos[0].x()), static_cast<int>(t.viewport_pos[1].x()), static_cast<int>(t.viewport_pos[2].x())});
    int max_x = std::max({static_cast<int>(t.viewport_pos[0].x()), static_cast<int>(t.viewport_pos[1].x()), static_cast<int>(t.viewport_pos[2].x())});
    int min_y = std::min({static_cast<int>(t.viewport_pos[0].y()), static_cast<int>(t.viewport_pos[1].y()), static_cast<int>(t.viewport_pos[2].y())});
    int max_y = std::max({static_cast<int>(t.viewport_pos[0].y()), static_cast<int>(t.viewport_pos[1].y()), static_cast<int>(t.viewport_pos[2].y())});

    // 遍历边界框内的每个像素
    for (int y = min_y; y <= max_y; ++y) {
        for (int x = min_x; x <= max_x; ++x) {
            // 检查像素 (x, y) 是否在三角形内部
            if (inside_triangle(x, y, t.viewport_pos)) {
                // 计算重心坐标
                auto [alpha, beta, gamma] = compute_barycentric_2d(float(x), float(y), t.viewport_pos);

                // 插值深度 (z 值)，使用投影校正算法
                float interpolated_reciprocal_w = alpha / t.world_pos[0].w() + beta / t.world_pos[1].w() + gamma / t.world_pos[2].w();
                float interpolated_depth = 1.0f / interpolated_reciprocal_w;

                // 插值其他顶点属性（如法线）
                Eigen::Vector3f interpolated_normal =
                    interpolate(alpha, beta, gamma, t.normal[0], t.normal[1], t.normal[2],
                                {t.world_pos[0].w(), t.world_pos[1].w(), t.world_pos[2].w()},
                                interpolated_depth);

                // 创建片段并推送到片段队列


                FragmentShaderPayload payload;
                payload.world_pos =
                    interpolate(alpha, beta, gamma,
                                Eigen::Vector3f(t.world_pos[0].x(), t.world_pos[0].y(),
                                                t.world_pos[0].z()), // 提取前三个分量
                                Eigen::Vector3f(t.world_pos[1].x(), t.world_pos[1].y(),
                                                t.world_pos[1].z()), // 提取前三个分量
                                Eigen::Vector3f(t.world_pos[2].x(), t.world_pos[2].y(),
                                                t.world_pos[2].z()), // 提取前三个分量
                                {t.world_pos[0].w(), t.world_pos[1].w(), t.world_pos[2].w()},
                                interpolated_depth);
                payload.world_normal = interpolated_normal.normalized();
                payload.x = x;
                payload.y = y;
                payload.depth = interpolated_depth;
                payload.color = Eigen::Vector3f(255.0f, 255.0f, 255.0f); // 假设为白色，可以根据实际情况修改

                // 上锁并推送片段
                std::unique_lock<std::mutex> lock(Context::rasterizer_queue_mutex);
                Context::rasterizer_output_queue.push(payload);
            }
        }
    }
}
