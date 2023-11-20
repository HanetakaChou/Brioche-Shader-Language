//
// Copyright (C) YuqiaoZhang(HanetakaChou)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//

#ifndef _BRX_BRDF_H_
#define _BRX_BRDF_H_ 1

// Avoid too small roughness to make the Trowbridge Reitz NDF work
// Real-Time Rendering Fourth Edition / 9.8.1 Normal Distribution Functions: "In the Disney principled shading model, Burley[214] exposes the roughness control to users as g = r2, where r is the user-interface roughness parameter value between 0 and 1."
// alpha = roughness * roughness;
// https://github.com/EpicGames/UnrealEngine/blob/4.27/Engine/Shaders/Private/CapsuleLightIntegrate.ush#L94
// NOTE: UE4 set "roughness" instead of "alpha" to 0.02, but alpha 0.02 is enough for us
static constexpr float const BRX_TROWBRIDGE_REITZ_ALPHA_MINIMUM = 0.02F;

// Prevent the NdotV to be zero
// https://github.com/EpicGames/UnrealEngine/blob/4.27/Engine/Shaders/Private/BRDF.ush#L34
static constexpr float const BRX_TROWBRIDGE_REITZ_NDOTV_MINIMUM = 1E-5F;

static constexpr float const BRX_BRDF_LENGTH_MINIMUM = 1E-5F;

static inline DirectX::XMFLOAT3 brx_trowbridge_reitz_sample_omega_h(DirectX::XMFLOAT2 const &xi, float const raw_alpha, DirectX::XMFLOAT3 const &raw_omega_o)
{
    // PBR Book V3: [Equation 8.12](https://pbr-book.org/3ed-2018/Reflection_Models/Microfacet_Models#MaskingandShadowing)
    // PBRT-V3: [TrowbridgeReitzDistribution::Sample_wh](https://github.com/mmp/pbrt-v3/blob/book/src/core/microfacet.cpp#L308)
    // PBR Book V4: [Equation 9.23](https://pbr-book.org/4ed/Reflection_Models/Roughness_Using_Microfacet_Theory#SamplingtheDistributionofVisibleNormals)
    // PBRT-V4: [TrowbridgeReitzDistribution::Sample_wm](https://github.com/mmp/pbrt-v4/blob/master/src/pbrt/util/scattering.h#L163)
    // UE4: [ImportanceSampleVisibleGGX](https://github.com/EpicGames/UnrealEngine/blob/4.27/Engine/Shaders/Private/MonteCarlo.ush#L380)
    // U3D: [SampleGGXVisibleNormal](https://github.com/Unity-Technologies/Graphics/blob/v10.8.1/com.unity.render-pipelines.core/ShaderLibrary/ImageBasedLighting.hlsl#L222)

    assert(raw_alpha >= BRX_TROWBRIDGE_REITZ_ALPHA_MINIMUM);
    float const alpha = std::max(BRX_TROWBRIDGE_REITZ_ALPHA_MINIMUM, raw_alpha);

    assert(raw_omega_o.z >= BRX_TROWBRIDGE_REITZ_NDOTV_MINIMUM);
    DirectX::XMFLOAT3 const omega_o(raw_omega_o.x, raw_omega_o.y, std::max(BRX_TROWBRIDGE_REITZ_NDOTV_MINIMUM, raw_omega_o.z));

    DirectX::XMFLOAT3 omega_o_hemisphere;
    {
        DirectX::XMFLOAT3 omega_o_hemisphere_non_unit(omega_o.x * alpha, omega_o.y * alpha, omega_o.z);
        DirectX::XMStoreFloat3(&omega_o_hemisphere, DirectX::XMVector3Normalize(DirectX::XMLoadFloat3(&omega_o_hemisphere_non_unit)));
    }

    DirectX::XMFLOAT3 T1;
    {
        DirectX::XMFLOAT3 T_1_raw = DirectX::XMFLOAT3(-omega_o_hemisphere.y, omega_o_hemisphere.x, 0.0F);
        float T_1_length_square = DirectX::XMVectorGetX(DirectX::XMVector3Dot(DirectX::XMLoadFloat3(&T_1_raw), DirectX::XMLoadFloat3(&T_1_raw)));
        if (T_1_length_square > BRX_BRDF_LENGTH_MINIMUM)
        {
            DirectX::XMStoreFloat3(&T1, DirectX::XMVectorScale(DirectX::XMLoadFloat3(&T_1_raw), 1.0F / std::sqrt(T_1_length_square)));
        }
        else
        {
            T1 = DirectX::XMFLOAT3(1.0F, 0.0F, 0.0F);
        }
    }

    DirectX::XMFLOAT3 T2;
    {
        DirectX::XMStoreFloat3(&T2, DirectX::XMVector3Cross(DirectX::XMLoadFloat3(&omega_o_hemisphere), DirectX::XMLoadFloat3(&T1)));
    }

    DirectX::XMFLOAT3 p;
    {
        float r = std::sqrt(std::max(0.0F, xi.x));
        float theta = 2.0F * DirectX::XM_PI * xi.y;

#if 0
        float cos_theta;
        float sin_theta;
        DirectX::XMScalarSinCos(&sin_theta,&cos_theta, theta);
#else
        float cos_theta = std::cos(theta);
        float sin_theta = std::sin(theta);
#endif

        float disk_x = r * cos_theta;
        float disk_y = r * sin_theta;

        float p_x = disk_x;

        float p_y;
        {
            float lerp_factor = (1.0F + omega_o_hemisphere.z) * 0.5F;
            p_y = std::sqrt(std::max(0.0F, 1.0F - disk_x * disk_x)) * (1.0F - lerp_factor) + disk_y * lerp_factor;
        }

        float p_z;
        {
            DirectX::XMFLOAT2 p_xy(p_x, p_y);
            p_z = std::sqrt(std::max(0.0F, 1.0F - DirectX::XMVectorGetX(DirectX::XMVector2Dot(DirectX::XMLoadFloat2(&p_xy), DirectX::XMLoadFloat2(&p_xy)))));
        }

        p = DirectX::XMFLOAT3(p_x, p_y, p_z);
    }

    DirectX::XMFLOAT3 n_h;
    {
        DirectX::XMStoreFloat3(&n_h, DirectX::XMVectorAdd(DirectX::XMVectorAdd(DirectX::XMVectorScale(DirectX::XMLoadFloat3(&T1), p.x), DirectX::XMVectorScale(DirectX::XMLoadFloat3(&T2), p.y)), DirectX::XMVectorScale(DirectX::XMLoadFloat3(&omega_o_hemisphere), p.z)));
    }

    DirectX::XMFLOAT3 omega_h;
    {
        assert(n_h.z >= BRX_TROWBRIDGE_REITZ_NDOTV_MINIMUM);
        DirectX::XMFLOAT3 omega_h_non_unit(n_h.x * alpha, n_h.y * alpha, std::max(BRX_TROWBRIDGE_REITZ_NDOTV_MINIMUM, n_h.z));
        DirectX::XMStoreFloat3(&omega_h, DirectX::XMVector3Normalize(DirectX::XMLoadFloat3(&omega_h_non_unit)));
    }
    return omega_h;
}

static inline float brx_trowbridge_reitz_pdf_omega_i(float const raw_alpha, float const raw_NdotV, float const NdotH)
{
    // VNDF = D * VdotH * G1 / NdotV
    //
    // PBR Book V3: [Equation 8.12](https://pbr-book.org/3ed-2018/Reflection_Models/Microfacet_Models#MaskingandShadowing)
    // PBRT-V3: [MicrofacetDistribution::Pdf](https://github.com/mmp/pbrt-v3/blob/book/src/core/microfacet.cpp#L339)
    // PBR Book V4: [Equation 9.23](https://pbr-book.org/4ed/Reflection_Models/Roughness_Using_Microfacet_Theory#SamplingtheDistributionofVisibleNormals)
    // PBRT-V4: [TrowbridgeReitzDistribution::PDF](https://github.com/mmp/pbrt-v4/blob/master/src/pbrt/util/scattering.h#L160)
    // UE4: [ImportanceSampleVisibleGGX](https://github.com/EpicGames/UnrealEngine/blob/4.27/Engine/Shaders/Private/MonteCarlo.ush#L380)
    // U3D: [SampleGGXVisibleNormal](https://github.com/Unity-Technologies/Graphics/blob/v10.8.1/com.unity.render-pipelines.core/ShaderLibrary/ImageBasedLighting.hlsl#L222)

    // PDF = VNDF / (4.0 * VdotH) = D * (G1 / NdotV) / 4.0
    //
    // PBR Book V3: [Figure 14.4](https://pbr-book.org/3ed-2018/Light_Transport_I_Surface_Reflection/Sampling_Reflection_Functions#MicrofacetBxDFs)
    // PBRT-V3: [MicrofacetReflection::Sample_f](https://github.com/mmp/pbrt-v3/blob/book/src/core/reflection.cpp#L413)
    // PBR Book V4: [Figure 9.30](https://pbr-book.org/4ed/Reflection_Models/Roughness_Using_Microfacet_Theory#x5-TheHalf-DirectionTransform)
    // PBRT-V4: [ConductorBxDF ::Sample_f](https://github.com/mmp/pbrt-v4/blob/master/src/pbrt/bxdfs.h#L317)

    float pdf;
    {
        assert(raw_alpha >= BRX_TROWBRIDGE_REITZ_ALPHA_MINIMUM);
        float const alpha = std::max(BRX_TROWBRIDGE_REITZ_ALPHA_MINIMUM, raw_alpha);

        assert(raw_NdotV >= BRX_TROWBRIDGE_REITZ_NDOTV_MINIMUM);
        float const NdotV = std::max(BRX_TROWBRIDGE_REITZ_NDOTV_MINIMUM, raw_NdotV);

        float alpha2 = alpha * alpha;

        float D;
        {
            // Real-Time Rendering Fourth Edition: Equation 9.41
            // PBR Book V3: [Figure 8.16](https://pbr-book.org/3ed-2018/Reflection_Models/Microfacet_Models#MicrofacetDistributionFunctions)
            // PBRT-V3: [TrowbridgeReitzDistribution::D](https://github.com/mmp/pbrt-v3/blob/book/src/core/microfacet.cpp#L153)
            // PBR Book V4: [Equation 9.16](https://pbr-book.org/4ed/Reflection_Models/Roughness_Using_Microfacet_Theory#TheMicrofacetDistribution)
            // PBRT-V4: [TrowbridgeReitzDistribution::D](https://github.com/mmp/pbrt-v4/blob/master/src/pbrt/util/scattering.h#L125)
            float term_h = NdotH * NdotH * (alpha2 - 1.0F) + 1.0F;
            D = alpha2 / (DirectX::XM_PI * term_h * term_h);
        }

        float G1_div_NdotV_div_4;
        {
            G1_div_NdotV_div_4 = 0.5F / (NdotV + std::sqrt(NdotV * (NdotV - NdotV * alpha2) + alpha2));
        }

        // float const PDF_MAXIMUM = (1.0F / (2.0F * DirectX::XM_PI * alpha2 * (NdotV + std::sqrt(NdotV * (NdotV - NdotV * alpha2) + alpha2))));
        // assert((D * G1_div_NdotV_div_4) <= PDF_MAXIMUM);
        // pdf = std::min(D * G1_div_NdotV_div_4, PDF_MAXIMUM);
        pdf = D * G1_div_NdotV_div_4;
    }

    return pdf;
}

static inline float brx_trowbridge_reitz_brdf_without_fresnel(float const raw_alpha, float const NdotH, float const raw_NdotV, float const NdotL)
{
    // https://pharr.org/matt/blog/2022/05/06/trowbridge-reitz

    assert(raw_alpha >= BRX_TROWBRIDGE_REITZ_ALPHA_MINIMUM);
    float const alpha = std::max(BRX_TROWBRIDGE_REITZ_ALPHA_MINIMUM, raw_alpha);

    assert(NdotH >= 0.0F);

    assert(raw_NdotV >= BRX_TROWBRIDGE_REITZ_NDOTV_MINIMUM);
    float const NdotV = std::max(BRX_TROWBRIDGE_REITZ_NDOTV_MINIMUM, raw_NdotV);

    assert(NdotL >= 0.0F);

    float alpha2 = alpha * alpha;

    float D;
    {
        // Real-Time Rendering Fourth Edition: Equation 9.41
        // PBR Book V3: [Figure 8.16](https://pbr-book.org/3ed-2018/Reflection_Models/Microfacet_Models#MicrofacetDistributionFunctions)
        // PBRT-V3: [TrowbridgeReitzDistribution::D](https://github.com/mmp/pbrt-v3/blob/book/src/core/microfacet.cpp#L153)
        // PBR Book V4: [Equation 9.16](https://pbr-book.org/4ed/Reflection_Models/Roughness_Using_Microfacet_Theory#TheMicrofacetDistribution)
        // PBRT-V4: [TrowbridgeReitzDistribution::D](https://github.com/mmp/pbrt-v4/blob/master/src/pbrt/util/scattering.h#L125)
        float term_h = NdotH * NdotH * (alpha2 - 1.0F) + 1.0F;
        D = alpha2 / (DirectX::XM_PI * term_h * term_h);
    }

    float V;
    {
        // Height-Correlated Trowbridge-Reitz

        // Lambda:
        // Equation 9.42 of Real-Time Rendering Fourth Edition
        // PBR Book V3: [Equation 8.13](https://pbr-book.org/3ed-2018/Reflection_Models/Microfacet_Models#MaskingandShadowing)
        // PBRT-V3: [TrowbridgeReitzDistribution::Lambda](https://github.com/mmp/pbrt-v3/blob/book/src/core/microfacet.cpp#L174)
        // PBR Book V4: [Equation9.22](https://pbr-book.org/4ed/Reflection_Models/Roughness_Using_Microfacet_Theory#TheMasking-ShadowingFunction)
        // PBRT-V4: [TrowbridgeReitzDistribution::Lambda](https://github.com/mmp/pbrt-v4/blob/master/src/pbrt/util/scattering.h#L143)
        // float lambda_v = 0.5 * (-1.0 + (1.0 / NdotV) * std::sqrt(alpha2 + (1.0 - alpha2) * NdotV * NdotV));
        // float lambda_l = 0.5 * (-1.0 + (1.0 / NdotL) * std::sqrt(alpha2 + (1.0 - alpha2) * NdotL * NdotL));

        // G2
        // Equation 9.31 of Real-Time Rendering Fourth Edition
        // PBR Book / 8.4.3 Masking and Shadowing: "A more accurate model can be derived assuming that microfacet visibility is more likely the higher up a given point on a microface"
        // G2 = 1.0 / (1.0 + Lambda(V) + Lambda(L)) = (2.0 * NoV * NoL) / (NoL * sqrt(alpha ^ 2 + (1.0 - alpha ^ 2) * NoV ^ 2) + NoV * sqrt(alpha ^ 2 + (1.0 - alpha ^ 2) * NoL ^ 2))
        // float G2 = 1.0 / (1.0 + lambda_v + lambda_l);

        // V = G2 / (4.0 * NoV * NoL) = 0.5 / (NoL * sqrt(alpha ^ 2 + (1.0 - alpha ^ 2) * NoV ^ 2) + NoV * sqrt(alpha ^ 2 + (1.0 - alpha ^ 2) * NoL ^ 2))
        float term_v = NdotL * std::sqrt(alpha2 + (1.0F - alpha2) * NdotV * NdotV);
        float term_l = NdotV * std::sqrt(alpha2 + (1.0F - alpha2) * NdotL * NdotL);

        V = (0.5F / (term_v + term_l));
    }

    // float const DV_MAXIMUM = (1.0F / (2.0F * DirectX::XM_PI * alpha2 * alpha * NdotV));
    // assert((D * V) <= DV_MAXIMUM);
    // float DV = std::min(D * V, DV_MAXIMUM);
    float DV = D * V;
    return DV;
}

static inline float brx_trowbridge_reitz_throughput_without_fresnel(float const raw_alpha, float const NdotL)
{
    assert(raw_alpha >= BRX_TROWBRIDGE_REITZ_ALPHA_MINIMUM);
    float const alpha = std::max(BRX_TROWBRIDGE_REITZ_ALPHA_MINIMUM, raw_alpha);

    assert(NdotL >= 0.0F);

    float alpha2 = alpha * alpha;

    assert((2.0F * NdotL / (NdotL + std::sqrt(NdotL * (NdotL - NdotL * alpha2) + alpha2))) <= 1.0F);
    float G2_div_G1 = std::min(2.0F * NdotL / (NdotL + std::sqrt(NdotL * (NdotL - NdotL * alpha2) + alpha2)), 1.0F);
    return G2_div_G1;
}

#endif
