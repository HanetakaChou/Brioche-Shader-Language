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

#ifndef _BRX_PACKED_VECTOR_H_
#define _BRX_PACKED_VECTOR_H_ 1

static inline uint32_t brx_FLOAT3_to_R15G15B2_SNORM(DirectX::XMFLOAT2 const &unpacked_input_xy, float unpacked_input_z)
{
    DirectX::XMFLOAT2 tmp_xy;
    float tmp_z;
    DirectX::XMStoreFloat2(&tmp_xy, DirectX::XMVectorScale(DirectX::XMVectorClamp(DirectX::XMLoadFloat2(&unpacked_input_xy), DirectX::XMVectorNegate(DirectX::XMVectorSplatOne()), DirectX::XMVectorSplatOne()), 16383.0F));
    tmp_z = std::copysign(1.0F, unpacked_input_z);

    int32_t element_x = (static_cast<int32_t>(tmp_xy.x) & 0x7fff);
    int32_t element_y = ((static_cast<int32_t>(tmp_xy.y) & 0x7fff) << 15);
    int32_t element_z = (static_cast<int32_t>(tmp_z) << 30);

    return static_cast<uint32_t>((element_x | element_y | element_z));
}

#endif
