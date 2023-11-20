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

#ifndef _BRX_LOW_DISCREPANCY_SEQUENCE_H_
#define _BRX_LOW_DISCREPANCY_SEQUENCE_H_ 1

// PBR Book V3: ["13.8.2 Quasi Monte Carlo"](https://www.pbr-book.org/3ed-2018/Monte_Carlo_Integration/Careful_Sample_Placement#QuasiMonteCarlo)
// PBR Book V4: ["8.2.2 Low Discrepancy and Quasi Monte Carlo"](https://pbr-book.org/4ed/Sampling_and_Reconstruction/Sampling_and_Integration#LowDiscrepancyandQuasiMonteCarlo)

static inline DirectX::XMFLOAT2 brx_hammersley_2d(uint32_t const sample_index, uint32_t const sample_count)
{
    // PBR Book V3: ["7.4.1 Hammersley and Halton Sequences"](https://www.pbr-book.org/3ed-2018/Sampling_and_Reconstruction/The_Halton_Sampler#HammersleyandHaltonSequences)
    // PBR Book V4: ["8.6.1 Hammersley and Halton Points"](https://pbr-book.org/4ed/Sampling_and_Reconstruction/Halton_Sampler#HammersleyandHaltonPoints)
    // UE: [Hammersley](https://github.com/EpicGames/UnrealEngine/blob/4.27/Engine/Shaders/Private/MonteCarlo.ush#L34)
    // U3D: [Hammersley2d](https://github.com/Unity-Technologies/Graphics/blob/v10.8.1/com.unity.render-pipelines.core/ShaderLibrary/Sampling/Hammersley.hlsl#L415)

    constexpr float const UINT32_MAX_PLUS_ONE = 4294967296.0F;

    float const xi_1 = static_cast<float>(sample_index) / static_cast<float>(sample_count);

#if defined(__GNUC__)
    // GCC or CLANG
    uint32_t const bit_field_reverse_sample_index = __builtin_bitreverse32(sample_index);
#elif defined(_MSC_VER)
#if defined(__clang__)
    // CLANG-CL
    uint32_t const bit_field_reverse_sample_index = __builtin_bitreverse32(sample_index);
#else
    // MSVC
    uint32_t bit_field_reverse_sample_index;
    {
        uint32_t bits = sample_index;
        bits = (bits << 16U) | (bits >> 16U);
        bits = ((bits & 0X00FF00FF) << 8U) | ((bits & 0XFF00FF00) >> 8U);
        bits = ((bits & 0X0F0F0F0F) << 4U) | ((bits & 0XF0F0F0F0) >> 4U);
        bits = ((bits & 0X33333333) << 2U) | ((bits & 0XCCCCCCCC) >> 2U);
        bits = ((bits & 0X55555555) << 1U) | ((bits & 0XAAAAAAAA) >> 1U);
        bit_field_reverse_sample_index = bits;
    }
#endif
#else
#error Unknown Compiler
#endif

    float const xi_2 = static_cast<float>(bit_field_reverse_sample_index) * (1.0F / static_cast<float>(UINT32_MAX_PLUS_ONE));

    return DirectX::XMFLOAT2(xi_1, xi_2);
}

#endif
