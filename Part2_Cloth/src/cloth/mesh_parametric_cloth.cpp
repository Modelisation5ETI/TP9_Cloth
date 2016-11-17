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


    //Gravity
    static vec3 const g (0.0f,0.0f,-9.81f);
    vec3 const g_normalized = g/N_total;
    for(int ku=0 ; ku<Nu ; ++ku)
    {
        for(int kv=0 ; kv<Nv ; ++kv)
        {
            force(ku,kv) = g_normalized;
        }
    }

    //*************************************************************//
    // TO DO, Calculer les forces s'appliquant sur chaque sommet
    //*************************************************************//

    //
    static float K = 20.0f;
    const float L_rest = 1.0f/Nu; //WARNING sqrt(2)/2 pour shearing et 2.0/Nu pour bending
    vec3 u;

    for(int ku=1 ; ku<Nu-1 ; ++ku)
    {
        for(int kv=1 ; kv<Nv-1 ; ++kv)
          {
          compute_neighbor_force(
            vec2(ku,kv),
            {vec2(ku+1,kv), vec2(ku - 1,kv), vec2(ku,kv + 1), vec2(ku,kv - 1)},
            K, L_rest);
          }
    }


    for(int kv=1 ; kv<Nv-1 ; ++kv)
      {
        compute_neighbor_force(vec2(0,kv),{vec2(1,kv),vec2(0,kv+1),vec2(0,kv-1)}, K, L_rest);


        compute_neighbor_force(vec2(Nu-1,kv),{vec2(Nu-2,kv),vec2(Nu-1,kv+1),vec2(Nu-1,kv-1)}, K, L_rest);

      }

    for(int ku=1 ; ku<Nu-1 ; ++ku)
      {
      compute_neighbor_force(vec2(ku,0),{vec2(ku,1),vec2(ku + 1,0),vec2(ku - 1,0)}, K, L_rest);


      compute_neighbor_force(vec2(ku,Nv-1),{vec2(ku,Nv-2),vec2(ku+1,Nv-1),vec2(ku-1,Nv-1)}, K, L_rest);


      }

    // Corners
    u = vertex(Nu-1, 0) - vertex(Nu-1, 1);
    force(Nu-1, 0) -= K * (norm(u)-L_rest) * u/norm(u);
    u = vertex(Nu-1, 0) - vertex(Nu-2, 0);
    force(Nu-1, 0) -= K * (norm(u)-L_rest) * u/norm(u);
//    force(Nu-1, 0) = force(Nu-1, 0)/2.0f;

    u = vertex(Nu-1,Nv-1) - vertex(Nu-2, Nv-1 );
    force(Nu-1,Nv-1) -= K * (norm(u)-L_rest) * u/norm(u);
    u = vertex(Nu-1,Nv-1) - vertex(Nu-1, Nv-2);
    force(Nu-1,Nv-1) -= K * (norm(u)-L_rest) * u/norm(u);
//    force(Nu-1,Nv-1) = force(Nu-1,Nv-1)/2.0f;

    u = vertex(0, 0) - vertex(1, 0);
    force(0, 0) -= K * (norm(u)-L_rest) * u/norm(u);
    u = vertex(0, 0) - vertex(0, 1);
    force(0, 0) -= K * (norm(u)-L_rest) * u/norm(u);
//    force(0, 0) = force(0, 0)/2.0f;

    u = vertex(0,Nv-1) - vertex(0, Nv-2);
    force(0,Nv-1) -= K * (norm(u)-L_rest) * u/norm(u);
    u = vertex(0,Nv-1) - vertex(1, Nv-1);
    force(0,Nv-1) -= K * (norm(u)-L_rest) * u/norm(u);
//    force(0,Nv-1) = force(0,Nv-1)/2.0f;



    // Point fixes
    force( 0, 0 ) = vec3( 0, 0, 0 );
    force( 0, Nv-1 ) = vec3( 0, 0, 0 );





    //*************************************************************//


}

void mesh_parametric_cloth::integration_step(float const dt)
{
    ASSERT_CPE(speed_data.size() == force_data.size(),"Incorrect size");
    ASSERT_CPE(static_cast<int>(speed_data.size()) == size_vertex(),"Incorrect size");


    int const Nu = size_u();
    int const Nv = size_v();
    //*************************************************************//
    // TO DO: Calculer l'integration numerique des positions au cours de l'intervalle de temps dt.
    //*************************************************************//
    static float const mu = 0.2f;
    for( int i = 0; i < Nu; i++ )
      {
      for( int j = 0; j < Nv; j++ )
        {
        speed(i,j) = (1-mu*dt)*speed(i,j)  + dt * force(i,j) ;
        vertex(i,j) = vertex(i,j) + dt * speed(i,j);
        }
      }


    //*************************************************************//


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
