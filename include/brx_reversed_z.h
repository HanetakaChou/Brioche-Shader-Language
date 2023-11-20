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

#ifndef _BRX_REVERSED_Z_H_
#define _BRX_REVERSED_Z_H_ 1

static inline DirectX::XMMATRIX XM_CALLCONV brx_DirectX_Math_Matrix_PerspectiveFovLH_ReversedZ(float FovAngleY, float AspectRatio, float NearZ, float FarZ) noexcept
{
    assert(NearZ > 0.f && FarZ > 0.f);
    assert(!DirectX::XMScalarNearEqual(FovAngleY, 0.0f, 0.00001f * 2.0f));
    assert(!DirectX::XMScalarNearEqual(AspectRatio, 0.0f, 0.00001f));
    assert(!DirectX::XMScalarNearEqual(FarZ, NearZ, 0.00001f));

    // [Reversed-Z](https://developer.nvidia.com/content/depth-precision-visualized)
    //
    // _  0  0  0
    // 0  _  0  0
    // 0  0  b  1
    // 0  0  a  0
    //
    // _  0  0  0
    // 0  _  0  0
    // 0  0 zb  z
    // 0  0  a
    //
    // z_prime = b + a / z
    //
    // Standard
    // 0 = b + a / nearz // z = nearz
    // 1 = b + a / farz // z = farz
    // a = farz * nearz / (nearz - farz)
    // b = - farz / (nearz - farz)
    //
    // Reversed-Z
    // 1 = b + a / nearz // z = nearz
    // 0 = b + a / farz // z = farz
    // a = farz * nearz / (farz - nearz)
    // b = - nearz / (farz - nearz)

    // DirectX::XMMatrixPerspectiveFovLH

    float SinFov;
    float CosFov;
    DirectX::XMScalarSinCos(&SinFov, &CosFov, 0.5F * FovAngleY);

    float Height = CosFov / SinFov;
    float Width = Height / AspectRatio;
    float b = -NearZ / (FarZ - NearZ);
    float a = (FarZ / (FarZ - NearZ)) * NearZ;

#if defined(_XM_ARM_NEON_INTRINSICS_)
    const float32x4_t Zero = vdupq_n_f32(0);

    DirectX::XMMATRIX M;
    M.r[0] = vsetq_lane_f32(Width, Zero, 0);
    M.r[1] = vsetq_lane_f32(Height, Zero, 1);
    M.r[2] = vsetq_lane_f32(b, DirectX::g_XMIdentityR3.v, 2);
    M.r[3] = vsetq_lane_f32(a, Zero, 2);
    return M;
#elif defined(_XM_SSE_INTRINSICS_)
    // __m128 _mm_shuffle_ps(__m128 lo,__m128 hi, _MM_SHUFFLE(hi3,hi2,lo1,lo0))
    // Interleave inputs into low 2 floats and high 2 floats of output. Basically
    // out[0] = lo[lo0];
    // out[1] = lo[lo1];
    // out[2] = hi[hi2];
    // out[3] = hi[hi3];

    // Note: This is recorded on the stack
    DirectX::XMVECTOR rMem = {
        Width,
        Height,
        b,
        a};

    // Copy from memory to SSE register
    DirectX::XMVECTOR vValues = rMem;
    DirectX::XMVECTOR vTemp = _mm_setzero_ps();
    // Copy x only
    vTemp = _mm_move_ss(vTemp, vValues);
    // CosFov / SinFov,0,0,0
    DirectX::XMMATRIX M;
    M.r[0] = vTemp;
    // 0,Height / AspectRatio,0,0
    vTemp = vValues;
    vTemp = _mm_and_ps(vTemp, DirectX::g_XMMaskY);
    M.r[1] = vTemp;
    // x=b,y=a,0,1.0f
    vTemp = _mm_setzero_ps();
    vValues = _mm_shuffle_ps(vValues, DirectX::g_XMIdentityR3, _MM_SHUFFLE(3, 2, 3, 2));
    // 0,0,b,1.0f
    vTemp = _mm_shuffle_ps(vTemp, vValues, _MM_SHUFFLE(3, 0, 0, 0));
    M.r[2] = vTemp;
    // 0,0,a,0.0f
    vTemp = _mm_shuffle_ps(vTemp, vValues, _MM_SHUFFLE(2, 1, 0, 0));
    M.r[3] = vTemp;
    return M;
#endif
}

static inline DirectX::XMMATRIX XM_CALLCONV brx_DirectX_Math_Matrix_PerspectiveFovRH_ReversedZ(float FovAngleY, float AspectRatio, float NearZ, float FarZ) noexcept
{
    assert(NearZ > 0.f && FarZ > 0.f);
    assert(!DirectX::XMScalarNearEqual(FovAngleY, 0.0f, 0.00001f * 2.0f));
    assert(!DirectX::XMScalarNearEqual(AspectRatio, 0.0f, 0.00001f));
    assert(!DirectX::XMScalarNearEqual(FarZ, NearZ, 0.00001f));

    // [Reversed-Z](https://developer.nvidia.com/content/depth-precision-visualized)
    //
    // _  0  0  0
    // 0  _  0  0
    // 0  0  b -1
    // 0  0  a  0
    //
    // _  0  0  0
    // 0  _  0  0
    // 0  0 zb  -z
    // 0  0  a
    //
    // z_prime = -b - a / z
    //
    // Standard
    // 0 = -b + a / nearz // z = -nearz
    // 1 = -b + a / farz // z = -farz
    // a = farz * nearz / (nearz - farz)
    // b = farz / (nearz - farz)
    //
    // Reversed-Z
    // 1 = -b + a / nearz // z = -nearz
    // 0 = -b + a / farz // z = -farz
    // a = farz * nearz / (farz - nearz)
    // b = nearz / (farz - nearz)

    // DirectX::XMMatrixPerspectiveFovRH

    float SinFov;
    float CosFov;
    DirectX::XMScalarSinCos(&SinFov, &CosFov, 0.5F * FovAngleY);

    float Height = CosFov / SinFov;
    float Width = Height / AspectRatio;
    float b = NearZ / (FarZ - NearZ);
    float a = (FarZ / (FarZ - NearZ)) * NearZ;

#if defined(_XM_ARM_NEON_INTRINSICS_)
    const float32x4_t Zero = vdupq_n_f32(0);

    DirectX::XMMATRIX M;
    M.r[0] = vsetq_lane_f32(Width, Zero, 0);
    M.r[1] = vsetq_lane_f32(Height, Zero, 1);
    M.r[2] = vsetq_lane_f32(b, DirectX::g_XMNegIdentityR3.v, 2);
    M.r[3] = vsetq_lane_f32(a, Zero, 2);
    return M;
#elif defined(_XM_SSE_INTRINSICS_)
    // __m128 _mm_shuffle_ps(__m128 lo,__m128 hi, _MM_SHUFFLE(hi3,hi2,lo1,lo0))
    // Interleave inputs into low 2 floats and high 2 floats of output. Basically
    // out[0] = lo[lo0];
    // out[1] = lo[lo1];
    // out[2] = hi[hi2];
    // out[3] = hi[hi3];

    // Note: This is recorded on the stack
    DirectX::XMVECTOR rMem = {
        Width,
        Height,
        b,
        a};

    // Copy from memory to SSE register
    DirectX::XMVECTOR vValues = rMem;
    DirectX::XMVECTOR vTemp = _mm_setzero_ps();
    // Copy x only
    vTemp = _mm_move_ss(vTemp, vValues);
    // CosFov / SinFov,0,0,0
    DirectX::XMMATRIX M;
    M.r[0] = vTemp;
    // 0,Height / AspectRatio,0,0
    vTemp = vValues;
    vTemp = _mm_and_ps(vTemp, DirectX::g_XMMaskY);
    M.r[1] = vTemp;
    // x=b,y=a,0,-1.0f
    vTemp = _mm_setzero_ps();
    vValues = _mm_shuffle_ps(vValues, DirectX::g_XMNegIdentityR3, _MM_SHUFFLE(3, 2, 3, 2));
    // 0,0,b,-1.0f
    vTemp = _mm_shuffle_ps(vTemp, vValues, _MM_SHUFFLE(3, 0, 0, 0));
    M.r[2] = vTemp;
    // 0,0,a,0.0f
    vTemp = _mm_shuffle_ps(vTemp, vValues, _MM_SHUFFLE(2, 1, 0, 0));
    M.r[3] = vTemp;
    return M;
#endif
}

#endif
