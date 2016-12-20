/*
**    TP CPE Lyon
**    Copyright (C) 2015 Damien Rohmer
**
**    This program is free software: you can redistribute it and/or modify
**    it under the terms of the GNU General Public License as published by
**    the Free Software Foundation, either version 3 of the License, or
**    (at your option) any later version.
**
**   This program is distributed in the hope that it will be useful,
**    but WITHOUT ANY WARRANTY; without even the implied warranty of
**    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**    GNU General Public License for more details.
**
**    You should have received a copy of the GNU General Public License
**    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "mesh_parametric_cloth.hpp"

#include "../lib/common/error_handling.hpp"
#include <cmath>

namespace cpe
{

void mesh_parametric_cloth::map_neighbor( const vec2& d1, const vec2& d2, const vec2& d3, const vec2& d4, const float& K, const float& L_rest )
{
  int const Nu = size_u();
  int const Nv = size_v();
  for(int ku=0 ; ku<Nu ; ++ku)
    {
    for(int kv=0 ; kv<Nv ; ++kv)
      {
      //Setup current neighborhood
      std::vector<vec2> neighbors;
      if( ku + d1.x() >= 0 && ku + d1.x() < Nu &&
          kv + d1.y() >= 0 && kv + d1.y() < Nv )
        {
        neighbors.push_back( vec2(ku,kv)+d1 );
        }
      if( ku + d2.x() >= 0 && ku + d2.x() < Nu &&
          kv + d2.y() >= 0 && kv + d2.y() < Nv )
        {
        neighbors.push_back( vec2(ku,kv)+d2 );
        }
      if( ku + d3.x() >= 0 && ku + d3.x() < Nu &&
          kv + d3.y() >= 0 && kv + d3.y() < Nv )
        {
        neighbors.push_back( vec2(ku,kv)+d3 );
        }
      if( ku + d4.x() >= 0 && ku + d4.x() < Nu &&
          kv + d4.y() >= 0 && kv + d4.y() < Nv )
        {
        neighbors.push_back( vec2(ku,kv)+d4 );
        }
      //Compute neighborhood forces
      compute_neighbor_force( vec2(ku,kv), neighbors, K, L_rest );
      }
    }
}

void mesh_parametric_cloth::compute_neighbor_force( const vec2& vertex_ij, std::vector<vec2> neighbors_ij, const float K, const float L_rest )
{
  for( unsigned int k=0; k < neighbors_ij.size(); k++  )
    {
    vec3 u = vertex(vertex_ij.x(),vertex_ij.y()) - vertex(neighbors_ij[k].x(),neighbors_ij[k].y());
    force(vertex_ij.x(),vertex_ij.y()) -= K * (norm(u)-L_rest) * u/norm(u);
    }
}

void mesh_parametric_cloth::update_force()
{
  int const Nu = size_u();
  int const Nv = size_v();
  int const N_total = Nu*Nv;
  ASSERT_CPE(static_cast<int>(force_data.size()) == Nu*Nv , "Error of size");

  // Gravity
  static vec3 const g (0.0f,0.0f,-9.81f);
  vec3 const g_normalized = g/N_total;

  // Wind
  vec3 windDirection( -1.5, 1, 0 );
  vec3 windDirection_normalized = windDirection/N_total;
  float windIntensity = 10.0f;

  // Apply forces
  for(int ku=0 ; ku<Nu ; ++ku)
    {
    for(int kv=0 ; kv<Nv ; ++kv)
      {
      force(ku,kv) = g_normalized;
      force(ku,kv) += windIntensity * dot( normal(ku,kv), windDirection_normalized ) * normal(ku,kv);
      }
    }

  //*************************************************************//
  // Compute spring neighbor forces
  //*************************************************************//
  // Structural
  map_neighbor( vec2(-1,0), vec2(0,-1), vec2(1,0), vec2(0,1), 30.0f, 1.0f/Nu );

  // Shear
  map_neighbor( vec2(-1,-1), vec2(1,-1), vec2(1,1), vec2(-1,1), 2.0f, sqrt(2)/(2*Nu) );

  // Bending
  map_neighbor( vec2(-2,0), vec2(0,-2), vec2(2,0), vec2(0,2), 25.0f, 2.0f/Nu );

  // Static Points
  force( 0, 0 ) = vec3( 0, 0, 0 );
  force( 0, Nv-1 ) = vec3( 0, 0, 0 );
}

void mesh_parametric_cloth::integration_step(float const dt)
{
  ASSERT_CPE(speed_data.size() == force_data.size(),"Incorrect size");
  ASSERT_CPE(static_cast<int>(speed_data.size()) == size_vertex(),"Incorrect size");

  int const Nu = size_u();
  int const Nv = size_v();
  //*************************************************************//
  // Numerical integration of positions and speeds along time
  //*************************************************************//
  static float const mu = 0.2f; //Damping coefficient
  for( int i = 0; i < Nu; i++ )
    {
    for( int j = 0; j < Nv; j++ )
      {
      //Euler explicit step
      speed(i,j) = (1-mu*dt)*speed(i,j)  + dt * force(i,j) ;
      vertex(i,j) = vertex(i,j) + dt * speed(i,j);

      //Collision with ground : WARNING: Hardcoded plane position
      if( vertex(i,j).z() < -1.101f +0.01f )
        {
        vertex(i,j).z() = -1.10f;
        speed(i,j).z() = -1.1f;
        }

      //Collision with sphere
      vec3 sphereCenter(0.4f,0.5f,-0.8f);
      float radius = 0.198f;
      vec3 normal = vertex(i,j)-sphereCenter;
      if( norm( normal ) < radius + 0.01 )
        {
        normal = normalized( normal );
        vertex(i,j) = sphereCenter + (radius + 0.01) * normal;
        //Set speed to 0 along normal
        speed(i,j) -= dot( normal, speed(i,j) ) * normal;
        }
      }
    }

  //security check (throw exception if divergence is detected)
  static float const LIMIT=30.0f;
  for(int ku=0 ; ku<Nu ; ++ku)
    {
    for(int kv=0 ; kv<Nv ; ++kv)
      {
      vec3 const& p = vertex(ku,kv);

      if( norm(p) > LIMIT )
        {
        throw exception_divergence("Divergence of the system",EXCEPTION_PARAMETERS_CPE);
        }
      }
    }

}

void mesh_parametric_cloth::set_plane_xy_unit(int const size_u_param,int const size_v_param)
{
    mesh_parametric::set_plane_xy_unit(size_u_param,size_v_param);

    int const N = size_u()*size_v();
    speed_data.resize(N);
    force_data.resize(N);
}

vec3 const& mesh_parametric_cloth::speed(int const ku,int const kv) const
{
    ASSERT_CPE(ku >= 0 , "Value ku ("+std::to_string(ku)+") should be >=0 ");
    ASSERT_CPE(ku < size_u() , "Value ku ("+std::to_string(ku)+") should be < size_u ("+std::to_string(size_u())+")");
    ASSERT_CPE(kv >= 0 , "Value kv ("+std::to_string(kv)+") should be >=0 ");
    ASSERT_CPE(kv < size_v() , "Value kv ("+std::to_string(kv)+") should be < size_v ("+std::to_string(size_v())+")");

    int const offset = ku + size_u()*kv;

    ASSERT_CPE(offset < static_cast<int>(speed_data.size()),"Internal error");

    return speed_data[offset];
}

vec3& mesh_parametric_cloth::speed(int const ku,int const kv)
{
    ASSERT_CPE(ku >= 0 , "Value ku ("+std::to_string(ku)+") should be >=0 ");
    ASSERT_CPE(ku < size_u() , "Value ku ("+std::to_string(ku)+") should be < size_u ("+std::to_string(size_u())+")");
    ASSERT_CPE(kv >= 0 , "Value kv ("+std::to_string(kv)+") should be >=0 ");
    ASSERT_CPE(kv < size_v() , "Value kv ("+std::to_string(kv)+") should be < size_v ("+std::to_string(size_v())+")");

    int const offset = ku + size_u()*kv;

    ASSERT_CPE(offset < static_cast<int>(speed_data.size()),"Internal error");

    return speed_data[offset];
}

vec3 const& mesh_parametric_cloth::force(int const ku,int const kv) const
{
    ASSERT_CPE(ku >= 0 , "Value ku ("+std::to_string(ku)+") should be >=0 ");
    ASSERT_CPE(ku < size_u() , "Value ku ("+std::to_string(ku)+") should be < size_u ("+std::to_string(size_u())+")");
    ASSERT_CPE(kv >= 0 , "Value kv ("+std::to_string(kv)+") should be >=0 ");
    ASSERT_CPE(kv < size_v() , "Value kv ("+std::to_string(kv)+") should be < size_v ("+std::to_string(size_v())+")");

    int const offset = ku + size_u()*kv;

    ASSERT_CPE(offset < static_cast<int>(force_data.size()),"Internal error");

    return force_data[offset];
}

vec3& mesh_parametric_cloth::force(int const ku,int const kv)
{
    ASSERT_CPE(ku >= 0 , "Value ku ("+std::to_string(ku)+") should be >=0 ");
    ASSERT_CPE(ku < size_u() , "Value ku ("+std::to_string(ku)+") should be < size_u ("+std::to_string(size_u())+")");
    ASSERT_CPE(kv >= 0 , "Value kv ("+std::to_string(kv)+") should be >=0 ");
    ASSERT_CPE(kv < size_v() , "Value kv ("+std::to_string(kv)+") should be < size_v ("+std::to_string(size_v())+")");

    int const offset = ku + size_u()*kv;

    ASSERT_CPE(offset < static_cast<int>(force_data.size()),"Internal error");

    return force_data[offset];
}


}
