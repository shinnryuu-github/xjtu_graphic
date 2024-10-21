#include "rasterizer_renderer.h"
#include "../utils/math.hpp"
#include <cstdio>
#include <iostream>

#ifdef _WIN32
#undef min
#undef max
#endif

using Eigen::Vector3f;
using Eigen::Vector4f;

// vertex shader
VertexShaderPayload vertex_shader(const VertexShaderPayload& payload)
{
    VertexShaderPayload output_payload = payload;

    // Vertex position transformation
    Eigen::Vector4f clip_position = Uniforms::MVP * payload.world_position;
    if (clip_position.w() != 0.0f) {
        clip_position /= clip_position.w();
    }
    output_payload.world_position = clip_position;
    // Viewport transformation
    float x = (clip_position.x() + 1.0f) * 0.5f * Uniforms::width;
    float y = (clip_position.y() + 1.0f) * 0.5f * Uniforms::height;
    output_payload.viewport_position = Eigen::Vector4f(x, y, clip_position.z(), 1.0f);

    // Vertex normal transformation
    Eigen::Vector4f normal_tmp = Eigen::Vector4f(payload.normal.x(), payload.normal.y(), payload.normal.z(), 0.0f); 
    normal_tmp = Uniforms::inv_trans_M * normal_tmp; 
    output_payload.normal = normal_tmp.head<3>().normalized();
    return output_payload;
}

Vector3f phong_fragment_shader(const FragmentShaderPayload& payload, const GL::Material& material,
                               const std::list<Light>& lights, const Camera& camera)
{
    // these lines below are just for compiling and can be deleted
    (void)payload;
    (void)material;
    (void)lights;
    (void)camera;
    // these lines above are just for compiling and can be deleted
    
    Vector3f result = {0, 0, 0};

    // 从材质中获取 ka, kd, ks
    // 设置环境光强度

    // 光照的贡献
    for (const auto& iter_light : lights) {
        // 计算光源方向
        Vector3f lightdir = (iter_light.position - payload.world_pos);
        lightdir.normalize();
        // 计算视点方向
        Vector3f viewdir = (camera.position - payload.world_pos);
        viewdir.normalize();
        // 计算半向量
        Vector3f half_vec = (lightdir + viewdir);
        half_vec.normalize();
        // 光的衰减
        float distance_squared = (iter_light.position - payload.world_pos).squaredNorm();
        float attenuation = iter_light.intensity / (distance_squared);  // 计算衰减

        // 环境光
        // 漫反射光
        float diff_factor = std::max(0.0f, (payload.world_normal).dot(lightdir));
        Vector3f diffuse = material.diffuse * attenuation * diff_factor; // 应用漫反射

        // 高光光
        Vector3f specular = material.specular * attenuation * std::pow(std::max(0.0f, (payload.world_normal).dot(half_vec)), material.shininess);

        std::cout << "distance_squared: " << distance_squared <<std::endl;
        std::cout << "light: " << attenuation <<std::endl;
        // std::cout << "payload.world_pos: " << payload.world_pos.transpose() <<std::endl;
        // std::cout << "iter_light.position: " << iter_light.position.transpose() <<std::endl;
        // 累加光照贡献
        result += material.ambient + diffuse + specular;
    }

    // 将渲染结果设置为最大阈值255
    result *= 255.f;
    // std::cout << "result: " << result.transpose() <<std::endl;
    return result;
}
