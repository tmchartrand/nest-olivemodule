/*
 *  hh_sdl_alpha_gap.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include "hh_sdl_alpha_gap.h"

// #ifdef HAVE_GSL

// C++ includes:
#include <cmath> // in case we need isnan() // fabs
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>

// Includes from libnestutil:
#include "numerics.h"

// Includes from nestkernel:
#include "exceptions.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"
#include "lockptrdatum.h"

using namespace nest;

nest::RecordablesMap< nest_tc::hh_sdl_alpha_gap >
  nest_tc::hh_sdl_alpha_gap::recordablesMap_;

namespace nest
{
// Override the create() method with one call to RecordablesMap::insert_()
// for each quantity to be recorded.
template <>
void
RecordablesMap< nest_tc::hh_sdl_alpha_gap >::create()
{
  // use standard names whereever you can for consistency!
  insert_( names::V_m,
    &nest_tc::hh_sdl_alpha_gap::get_y_elem_< nest_tc::hh_sdl_alpha_gap::State_::V_M > );
  insert_( names::I_syn_ex,
    &nest_tc::hh_sdl_alpha_gap::get_y_elem_< nest_tc::hh_sdl_alpha_gap::State_::I_EXC > );
  insert_( names::I_syn_in,
    &nest_tc::hh_sdl_alpha_gap::get_y_elem_< nest_tc::hh_sdl_alpha_gap::State_::I_INH > );
  // insert_( names::Act_m,
  //   &nest_tc::hh_sdl_alpha_gap::get_y_elem_< nest_tc::hh_sdl_alpha_gap::State_::M > );
  insert_( names::Act_h,
    &nest_tc::hh_sdl_alpha_gap::get_y_elem_< nest_tc::hh_sdl_alpha_gap::State_::H > );
  insert_( names::Inact_n,
    &nest_tc::hh_sdl_alpha_gap::get_y_elem_< nest_tc::hh_sdl_alpha_gap::State_::N > );
  // insert_( names::Inact_p,
  //   &nest_tc::hh_sdl_alpha_gap::get_y_elem_< nest_tc::hh_sdl_alpha_gap::State_::P > );
  // TC
  // insert_( names::l,
  //   &hh_sdl_alpha_gap::get_y_elem_< hh_sdl_alpha_gap::State_::L > );
  // insert_( names::q,
  //   &hh_sdl_alpha_gap::get_y_elem_< hh_sdl_alpha_gap::State_::Q > );
  // insert_( names::r,
  //   &hh_sdl_alpha_gap::get_y_elem_< hh_sdl_alpha_gap::State_::R > );
  // insert_( names::s,
  //   &hh_sdl_alpha_gap::get_y_elem_< hh_sdl_alpha_gap::State_::S > );

}
} // namespace nest

extern "C" int
nest_tc::hh_sdl_alpha_gap_dynamics( double time,
  const double y[],
  double f[],
  void* pnode )
{
  // a shorthand
  typedef nest_tc::hh_sdl_alpha_gap::State_ S;

  // get access to node so we can almost work as in a member function
  assert( pnode );
  const nest_tc::hh_sdl_alpha_gap& node =
    *( reinterpret_cast< nest_tc::hh_sdl_alpha_gap* >( pnode ) );

  // y[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.y[].

  // The following code is verbose for the sake of clarity. We assume that a
  // good compiler will optimize the verbosity away ...

  // shorthand for state variables
  const double& V = y[ S::V_M ];
  const double& Vd = y[ S::V_D ];
  const double& k = y[ S::K ];
  const double& l = y[ S::L ];
  const double& h = y[ S::H ];
  const double& n = y[ S::N ];
  const double& q = y[ S::Q ];
  const double& r = y[ S::R ];
  const double& s = y[ S::S ];
  const double& cca = y[ S::CCA ];
  const double& dI_ex = y[ S::DI_EXC ];
  const double& I_ex = y[ S::I_EXC ];
  const double& dI_in = y[ S::DI_INH ];
  const double& I_in = y[ S::I_INH ];

  // shorthand for params
  const double& E_Ca = node.P_.E_Ca;
  const double& E_Na = node.P_.E_Na;    
  const double& E_K = node.P_.E_K;      
  const double& E_H = node.P_.E_H;
  const double& E_L = node.P_.E_L;      
  const double& g_L = node.P_.g_L;      
  const double& g_Cal = node.P_.g_Cal;
  const double& g_Na = node.P_.g_Na;  
  const double& g_K = node.P_.g_K;
  const double& g_H = node.P_.g_H;
  const double& g_Int = node.P_.g_Int;  
  const double& g_Cah = node.P_.g_Cah;
  const double& g_Ld = node.P_.g_Ld;
  const double& g_Kca = node.P_.g_Kca;  
  const double& p_Ds = node.P_.p_Ds;

  
//Somatic

//Na
const double a_m =  0.1*(V+41)/(1-std::exp(-(V+41)/10));
const double b_m = 9.0*std::exp(-(V+66)/20);
const double mlim =  a_m/(a_m+b_m);
const double a_h =  5.0*std::exp(-(V+60)/15);
const double b_h =  (V+50)/(1-std::exp(-(V+50)/10));
// b_h =  (V+48)/(1-std::exp(-(V+48)/3)); 2004
const double t_h =  170.0/(a_h + b_h);
const double hlim =  a_h/(a_h+b_h);
const double I_na = g_Na*mlim*mlim*mlim*h*(V-E_Na);

//K_dr
const double a_n = (V+41)/(1-std::exp(-(V+41)/10));
const double b_n = 12.5*std::exp(-(V+51)/80);
const double nlim = a_n/(a_n+b_n);
const double t_n = 5.0/(a_n+b_n);
const double I_k = g_K*n*n*n*n*(V-E_K);

//Ca_l
const double klim = 1.0/(1+std::exp(-(V+61)/4.2));
const double t_k = 5;//1/5 dG? (=0 for manor)
const double llim = 1.0/(1+std::exp((V+85.5)/8.5));//
const double t_l = 35 + 20*std::exp((V+160)/30)/(1+std::exp((V+84)/7.3));
// t_l = 40 + 30*std::exp((V+160)/30)/(1+std::exp((V+84)/8.3));// manor
const double I_cal = g_Cal*k*k*k*l*(V-E_Ca);

//H
const double qlim = 1.0/(1+std::exp((V+75)/5.5));
// qlim = 1.0/(1+std::exp((V+80)/4));//vdg
const double t_q = 1.0/(std::exp(-0.086*V-14.6) + std::exp(0.07*V-1.87));
const double I_h = g_H*q*(V - E_H);

//leak, exchange
const double I_l = g_L*(V-E_L);
const double I_ds = g_Int/p_Ds*(V-Vd);

//Dendritic

//leak, exchange
const double I_sd = g_Int/(1-p_Ds)*(Vd-V);
const double I_ld = g_Ld*(Vd - E_L);

//Ca_h
const double a_r = 1.6/(1+std::exp(-(Vd-5)/14));
// a_r = 1.7/(1+std::exp(-(Vd+5)/13.9));//from deGruijl
const double b_r = -0.02*(Vd+8.5)/(1-std::exp((Vd+8.5)/5)); //*** must be -
const double rlim = a_r/(a_r+b_r);
const double t_r = 1.0/(a_r+b_r);//5/1? (dG)
const double I_cah = g_Cah*r*r*(Vd-E_Ca);

//K_Ca
const double a_s = std::min(0.00002*cca,0.01);
const double b_s = 0.015;
const double t_s = 1.0/(a_s+b_s);
const double slim = a_s/(a_s+b_s);
const double I_kca = g_Kca*s*(Vd-E_K);



  // set I_gap depending on interpolation order
  double gap = 0.0;

  const double t = time / node.B_.step_;

  switch ( kernel().simulation_manager.get_wfr_interpolation_order() )
  {
  case 0:
    gap = -node.B_.sumj_g_ij_ * V
      + node.B_.interpolation_coefficients[ node.B_.lag_ ];
    break;

  case 1:
    gap = -node.B_.sumj_g_ij_ * V
      + node.B_.interpolation_coefficients[ node.B_.lag_ * 2 + 0 ]
      + node.B_.interpolation_coefficients[ node.B_.lag_ * 2 + 1 ] * t;
    break;

  case 3:
    gap = -node.B_.sumj_g_ij_ * V
      + node.B_.interpolation_coefficients[ node.B_.lag_ * 4 + 0 ]
      + node.B_.interpolation_coefficients[ node.B_.lag_ * 4 + 1 ] * t
      + node.B_.interpolation_coefficients[ node.B_.lag_ * 4 + 2 ] * t * t
      + node.B_.interpolation_coefficients[ node.B_.lag_ * 4 + 3 ] * t * t * t;
    break;

  default:
    throw BadProperty( "Interpolation order must be 0, 1, or 3." );
  }

  const double I_gap = gap;

  // V dot -- synaptic input are currents, inhib current is negative
  f[ S::V_M ] = ( -(I_l+I_cal+I_na+I_k+I_h+I_ds)
    + node.B_.I_stim_ + node.P_.I_e ) / node.P_.C_m;

  f[ S::V_D ] = ( -(I_sd+I_ld+I_cah+I_kca)
     + I_ex + I_in + I_gap ) / node.P_.C_m;

  
  f[ S::K ] = -(k - klim)/t_k;
  f[ S::L ] = -(l - llim)/t_l;
  f[ S::H ] = -(h - hlim)/t_h;
  f[ S::N ] = -(n - nlim)/t_n;
  f[ S::Q ] = -(q - qlim)/t_q;
  f[ S::R ] = -(r - rlim)/t_r;
  f[ S::S ] = -(s - slim)/t_s;

  f[S::CCA] = -3.0*I_cah - 0.075*cca;

  // synapses: alpha functions
  f[ S::DI_EXC ] = -dI_ex / node.P_.tau_synE;
  f[ S::I_EXC ] = dI_ex - ( I_ex / node.P_.tau_synE );
  f[ S::DI_INH ] = -dI_in / node.P_.tau_synI;
  f[ S::I_INH ] = dI_in - ( I_in / node.P_.tau_synI );

  return GSL_SUCCESS;
}



/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

nest_tc::hh_sdl_alpha_gap::Parameters_::Parameters_()
  : tau_synE( 0.2 ) // ms
  , tau_synI( 2.0 ) // ms
  , I_e( 0.0 )      // pA 
  , C_m( 1. )     
  , t_ref_( 5. )   // ms

  , g_L( 0.015 )   
  , g_Ld( 0.015 )   
  , g_Na( 70. )   
  , g_K(18.)
  , g_Cal( 1.5)
  , g_H(1.5)
  , g_Cah(4.)
  , g_Kca( 35.)
  , g_Int( 0.13 )
  , p_Ds( 0.2 )

  , E_L( -63. )     // mV
  , E_Ca(120.)
  , E_Na( 55. )    // mV
  , E_K( -75. )    // mV
  , E_H(-43.)
{
}

nest_tc::hh_sdl_alpha_gap::State_::State_( const Parameters_& )
  : r_( 0 )
{
  y_[ 0 ] = -57.;
  y_[ 1 ] = -57.;
  y_[ 2 ] = 5.;
  for ( size_t i = 3; i < STATE_VEC_SIZE; ++i )
  {
    y_[ i ] = 0;
  }

  // equilibrium values for (in)activation variables
  // todo
}

nest_tc::hh_sdl_alpha_gap::State_::State_( const State_& s )
  : r_( s.r_ )
{
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
  {
    y_[ i ] = s.y_[ i ];
  }
}

nest_tc::hh_sdl_alpha_gap::State_& nest_tc::hh_sdl_alpha_gap::State_::operator=(
  const State_& s )
{
  assert( this != &s ); // would be bad logical error in program
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
  {
    y_[ i ] = s.y_[ i ];
  }
  r_ = s.r_;
  return *this;
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
nest_tc::hh_sdl_alpha_gap::Parameters_::get( DictionaryDatum& d ) const
{
  def< double >( d, names::t_ref, t_ref_ );
  def< double >( d, names::g_Na, g_Na );
  // def< double >( d, names::g_Kv1, g_Kv1 );
  // def< double >( d, names::g_Kv3, g_Kv3 );
  def< double >( d, names::g_L, g_L );
  def< double >( d, names::E_Na, E_Na );
  def< double >( d, names::E_K, E_K );
  def< double >( d, names::E_L, E_L );
  def< double >( d, names::C_m, C_m );
  def< double >( d, names::tau_syn_ex, tau_synE );
  def< double >( d, names::tau_syn_in, tau_synI );
  def< double >( d, names::I_e, I_e );
}

void
nest_tc::hh_sdl_alpha_gap::Parameters_::set( const DictionaryDatum& d )
{
  updateValue< double >( d, names::t_ref, t_ref_ );
  updateValue< double >( d, names::C_m, C_m );
  updateValue< double >( d, names::g_Na, g_Na );
  updateValue< double >( d, names::E_Na, E_Na );
  // updateValue< double >( d, names::g_Kv1, g_Kv1 );
  // updateValue< double >( d, names::g_Kv3, g_Kv3 );
  updateValue< double >( d, names::E_K, E_K );
  updateValue< double >( d, names::g_L, g_L );
  updateValue< double >( d, names::E_L, E_L );

  updateValue< double >( d, names::tau_syn_ex, tau_synE );
  updateValue< double >( d, names::tau_syn_in, tau_synI );

  updateValue< double >( d, names::I_e, I_e );
  if ( C_m <= 0 )
  {
    throw BadProperty( "Capacitance must be strictly positive." );
  }
  if ( t_ref_ < 0 )
  {
    throw BadProperty( "Refractory time cannot be negative." );
  }
  if ( tau_synE <= 0 || tau_synI <= 0 )
  {
    throw BadProperty( "All time constants must be strictly positive." );
  }
  if ( g_Na < 0 || g_L < 0 )
  {
    throw BadProperty( "All conductances must be non-negative." );
  }
}

void
nest_tc::hh_sdl_alpha_gap::State_::get( DictionaryDatum& d ) const
{
  def< double >( d, names::V_m, y_[ V_M ] );
  // def< double >( d, names::Act_m, y_[ M ] );
  def< double >( d, names::Act_h, y_[ H ] );
  // def< double >( d, names::Inact_n, y_[ N ] );
  // def< double >( d, names::Inact_p, y_[ P ] );
}

void
nest_tc::hh_sdl_alpha_gap::State_::set( const DictionaryDatum& d )
{
  updateValue< double >( d, names::V_m, y_[ V_M ] );
  // updateValue< double >( d, names::Act_m, y_[ M ] );
  updateValue< double >( d, names::Act_h, y_[ H ] );
  // updateValue< double >( d, names::Inact_n, y_[ N ] );
  // updateValue< double >( d, names::Inact_p, y_[ P ] );
  // if ( y_[ M ] < 0 || y_[ H ] < 0 || y_[ N ] < 0 || y_[ P ] < 0 )
  // {
  //   throw BadProperty( "All (in)activation variables must be non-negative." );
  // }
}

nest_tc::hh_sdl_alpha_gap::Buffers_::Buffers_( hh_sdl_alpha_gap& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

nest_tc::hh_sdl_alpha_gap::Buffers_::Buffers_( const Buffers_&,
  hh_sdl_alpha_gap& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

nest_tc::hh_sdl_alpha_gap::hh_sdl_alpha_gap()
  : Archiving_Node()
  , P_()
  , S_( P_ )
  , B_( *this )
{
  recordablesMap_.create();
  Node::set_node_uses_wfr( kernel().simulation_manager.use_wfr() );
}

nest_tc::hh_sdl_alpha_gap::hh_sdl_alpha_gap( const hh_sdl_alpha_gap& n )
  : Archiving_Node( n )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
  Node::set_node_uses_wfr( kernel().simulation_manager.use_wfr() );
}

nest_tc::hh_sdl_alpha_gap::~hh_sdl_alpha_gap()
{
  // GSL structs may not have been allocated, so we need to protect destruction
  if ( B_.s_ )
  {
    gsl_odeiv_step_free( B_.s_ );
  }
  if ( B_.c_ )
  {
    gsl_odeiv_control_free( B_.c_ );
  }
  if ( B_.e_ )
  {
    gsl_odeiv_evolve_free( B_.e_ );
  }
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
nest_tc::hh_sdl_alpha_gap::init_state_( const Node& proto )
{
  const hh_sdl_alpha_gap& pr = downcast< hh_sdl_alpha_gap >( proto );
  S_ = pr.S_;
}

void
nest_tc::hh_sdl_alpha_gap::init_buffers_()
{
  B_.spike_exc_.clear(); // includes resize
  B_.spike_inh_.clear(); // includes resize
  B_.currents_.clear();  // includes resize

  // allocate strucure for gap events here
  // function is called from Scheduler::prepare_nodes() before the
  // first call to update
  // so we already know which interpolation scheme to use according
  // to the properties of this neurons
  // determine size of structure depending on interpolation scheme
  // and unsigned int Scheduler::min_delay() (number of simulation time steps
  // per min_delay step)

  // resize interpolation_coefficients depending on interpolation order
  const size_t quantity = kernel().connection_manager.get_min_delay()
    * ( kernel().simulation_manager.get_wfr_interpolation_order() + 1 );

  B_.interpolation_coefficients.resize( quantity, 0.0 );

  B_.last_y_values.resize( kernel().connection_manager.get_min_delay(), 0.0 );

  B_.sumj_g_ij_ = 0.0;

  Archiving_Node::clear_history();

  B_.logger_.reset();

  B_.step_ = Time::get_resolution().get_ms();
  B_.IntegrationStep_ = B_.step_;

  if ( B_.s_ == 0 )
  {
    B_.s_ =
      gsl_odeiv_step_alloc( gsl_odeiv_step_rkf45, State_::STATE_VEC_SIZE );
  }
  else
  {
    gsl_odeiv_step_reset( B_.s_ );
  }

  if ( B_.c_ == 0 )
  {
    B_.c_ = gsl_odeiv_control_y_new( 1e-6, 0.0 );
  }
  else
  {
    gsl_odeiv_control_init( B_.c_, 1e-6, 0.0, 1.0, 0.0 );
  }

  if ( B_.e_ == 0 )
  {
    B_.e_ = gsl_odeiv_evolve_alloc( State_::STATE_VEC_SIZE );
  }
  else
  {
    gsl_odeiv_evolve_reset( B_.e_ );
  }

  B_.sys_.function = hh_sdl_alpha_gap_dynamics;
  B_.sys_.jacobian = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params = reinterpret_cast< void* >( this );

  B_.I_stim_ = 0.0;
}

void
nest_tc::hh_sdl_alpha_gap::calibrate()
{
  // ensures initialization in case mm connected after Simulate
  B_.logger_.init();

  V_.PSCurrInit_E_ = 1.0 * numerics::e / P_.tau_synE;
  V_.PSCurrInit_I_ = 1.0 * numerics::e / P_.tau_synI;
  V_.RefractoryCounts_ = Time( Time::ms( P_.t_ref_ ) ).get_steps();
  // since t_ref_ >= 0, this can only fail in error
  assert( V_.RefractoryCounts_ >= 0 );
}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

bool
nest_tc::hh_sdl_alpha_gap::update_( Time const& origin,
  const long from,
  const long to,
  const bool wfr_update )
{

  assert(
    to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
  assert( from < to );

  bool done = true;
  const size_t interpolation_order =
    kernel().simulation_manager.get_wfr_interpolation_order();
  const double wfr_tol = kernel().simulation_manager.get_wfr_tol();

  // allocate memory to store the new interpolation coefficients
  // to be sent by gap event
  const size_t quantity =
    kernel().connection_manager.get_min_delay() * ( interpolation_order + 1 );
  std::vector< double > new_coefficients( quantity, 0.0 );

  // parameters needed for piecewise interpolation
  double y_i = 0.0, y_ip1 = 0.0, hf_i = 0.0, hf_ip1 = 0.0;
  double f_temp[ State_::STATE_VEC_SIZE ];

  for ( long lag = from; lag < to; ++lag )
  {

    // B_.lag is needed by hh_sdl_alpha_gap_dynamics to
    // determine the current section
    B_.lag_ = lag;

    if ( wfr_update )
    {
      y_i = S_.y_[ State_::V_M ];
      if ( interpolation_order == 3 )
      {
        hh_sdl_alpha_gap_dynamics(
          0, S_.y_, f_temp, reinterpret_cast< void* >( this ) );
        hf_i = B_.step_ * f_temp[ State_::V_M ];
      }
    }

    double t = 0.0;
    const double U_old = S_.y_[ State_::V_M ];

    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t; this is of advantage
    // for a consistent and efficient integration across subsequent
    // simulation intervals
    while ( t < B_.step_ )
    {
      const int status = gsl_odeiv_evolve_apply( B_.e_,
        B_.c_,
        B_.s_,
        &B_.sys_,             // system of ODE
        &t,                   // from t
        B_.step_,             // to t <= step
        &B_.IntegrationStep_, // integration step size
        S_.y_ );              // neuronal state
      if ( status != GSL_SUCCESS )
      {
        throw GSLSolverFailure( get_name(), status );
      }
    }

    if ( not wfr_update )
    {
      S_.y_[ State_::DI_EXC ] +=
        B_.spike_exc_.get_value( lag ) * V_.PSCurrInit_E_;
      S_.y_[ State_::DI_INH ] +=
        B_.spike_inh_.get_value( lag ) * V_.PSCurrInit_I_;
      // sending spikes: crossing 0 mV, pseudo-refractoriness and local
      // maximum...
      // refractory?
      if ( S_.r_ > 0 )
      {
        --S_.r_;
      }
      else
        // (    threshold    &&     maximum       )
        if ( S_.y_[ State_::V_M ] >= 0 && U_old > S_.y_[ State_::V_M ] )
      {
        S_.r_ = V_.RefractoryCounts_;

        set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );

        SpikeEvent se;
        kernel().event_delivery_manager.send( *this, se, lag );
      }

      // log state data
      B_.logger_.record_data( origin.get_steps() + lag );

      // set new input current
      B_.I_stim_ = B_.currents_.get_value( lag );
    }
    else // if(wfr_update)
    {
      S_.y_[ State_::DI_EXC ] +=
        B_.spike_exc_.get_value_wfr_update( lag ) * V_.PSCurrInit_E_;
      S_.y_[ State_::DI_INH ] +=
        B_.spike_inh_.get_value_wfr_update( lag ) * V_.PSCurrInit_I_;
      // check deviation from last iteration
      done = ( fabs( S_.y_[ State_::V_M ] - B_.last_y_values[ lag ] )
               <= wfr_tol ) && done;
      B_.last_y_values[ lag ] = S_.y_[ State_::V_M ];

      // update different interpolations

      // constant term is the same for each interpolation order
      new_coefficients[ lag * ( interpolation_order + 1 ) + 0 ] = y_i;

      switch ( interpolation_order )
      {
      case 0:
        break;

      case 1:
        y_ip1 = S_.y_[ State_::V_M ];

        new_coefficients[ lag * ( interpolation_order + 1 ) + 1 ] = y_ip1 - y_i;
        break;

      case 3:
        y_ip1 = S_.y_[ State_::V_M ];
        hh_sdl_alpha_gap_dynamics(
          B_.step_, S_.y_, f_temp, reinterpret_cast< void* >( this ) );
        hf_ip1 = B_.step_ * f_temp[ State_::V_M ];

        new_coefficients[ lag * ( interpolation_order + 1 ) + 1 ] = hf_i;
        new_coefficients[ lag * ( interpolation_order + 1 ) + 2 ] =
          -3 * y_i + 3 * y_ip1 - 2 * hf_i - hf_ip1;
        new_coefficients[ lag * ( interpolation_order + 1 ) + 3 ] =
          2 * y_i - 2 * y_ip1 + hf_i + hf_ip1;
        break;

      default:
        throw BadProperty( "Interpolation order must be 0, 1, or 3." );
      }
    }


  } // end for-loop

  // if not wfr_update perform constant extrapolation and reset last_y_values
  if ( not wfr_update )
  {
    for ( long temp = from; temp < to; ++temp )
    {
      new_coefficients[ temp * ( interpolation_order + 1 ) + 0 ] =
        S_.y_[ State_::V_M ];
    }

    B_.last_y_values.clear();
    B_.last_y_values.resize( kernel().connection_manager.get_min_delay(), 0.0 );
  }

  // Send gap-event
  GapJunctionEvent ge;
  ge.set_coeffarray( new_coefficients );
  kernel().event_delivery_manager.send_secondary( *this, ge );

  // Reset variables
  B_.sumj_g_ij_ = 0.0;
  B_.interpolation_coefficients.clear();
  B_.interpolation_coefficients.resize( quantity, 0.0 );

  return done;
}

void
nest_tc::hh_sdl_alpha_gap::handle( SpikeEvent& e )
{
  assert( e.get_delay() > 0 );

  if ( e.get_weight() > 0.0 )
  {
    B_.spike_exc_.add_value( e.get_rel_delivery_steps(
                               kernel().simulation_manager.get_slice_origin() ),
      e.get_weight() * e.get_multiplicity() );
  }
  else
  {
    B_.spike_inh_.add_value( e.get_rel_delivery_steps(
                               kernel().simulation_manager.get_slice_origin() ),
      e.get_weight() * e.get_multiplicity() );
  } // current input, keep negative weight
}

void
nest_tc::hh_sdl_alpha_gap::handle( CurrentEvent& e )
{
  assert( e.get_delay() > 0 );

  const double c = e.get_current();
  const double w = e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(
    e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ),
    w * c );
}

void
nest_tc::hh_sdl_alpha_gap::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}

void
nest_tc::hh_sdl_alpha_gap::handle( GapJunctionEvent& e )
{

  B_.sumj_g_ij_ += e.get_weight();

  size_t i = 0;
  std::vector< unsigned int >::iterator it = e.begin();
  // The call to get_coeffvalue( it ) in this loop also advances the iterator it
  while ( it != e.end() )
  {
    B_.interpolation_coefficients[ i ] +=
      e.get_weight() * e.get_coeffvalue( it );
    i++;
  }
}

// #endif // HAVE_GSL
